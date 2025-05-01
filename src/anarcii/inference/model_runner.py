import torch

from anarcii.input_data_processing import TokenisedSequence
from anarcii.input_data_processing.tokeniser import NumberingTokeniser

from .model_loader import Loader
from .utils import build_inward_list, dataloader

# NEED TO COME BACK TO THIS CODE AND LOOK AT THE TRY EXCEPT LOOPS....
# SOMETHING SHOULD BE MODIFIED TO REDUCE THEM....


# A cutoff score to consider a sequence as well numbered by the language model.
CUTOFF_SCORE = 15


class ModelRunner:
    """
    This class orchestrates the auto-regressive inference steps.
    It takes the list of tokensised seqs, pads then batches into a dataloader.
    The dataloader is iterated through and batches processed.

    Numbering is predicted in parallel for each batch in an autoregressive inference
    loop that iterates for as many steps as the longest seq per batch.

    Predicted tokens are then tranlated back to corresponding number and instertion
    labels.

    Each sequence is processed one by one in a for loop. In the for loop the chain call
    and sequence score are identifed, then based on this score and the number of
    integer labels (non insertion [x] labels) then the full sequence is processed.

    Processin invovles working out where seq numbering begins (detect [<skip>] tokens),
    if some numbering is missed from the begininng or end (backfill and forward fill)
    and when we have insertions [X] translating them to a number based on the labels
    either side [X] [X] [X] >> (111, A) (112, B) (112, A).

    Finally the sequence is appended with empty labels to ensure that the returned
    result contains all integer labels from 1 to 128.

    """

    def __init__(self, sequence_type, mode, batch_size, device, verbose):
        self.type = sequence_type.lower()
        self.mode = mode.lower()
        self.batch_size = batch_size
        self.device = device
        self.verbose = verbose

        if self.type in ["antibody", "shark"]:
            self.sequence_tokeniser = NumberingTokeniser("protein_antibody")
            self.number_tokeniser = NumberingTokeniser("number_antibody")

        elif self.type == "tcr":
            self.sequence_tokeniser = NumberingTokeniser("protein_tcr")
            self.number_tokeniser = NumberingTokeniser("number_tcr")
        else:
            raise ValueError(f"Invalid model type: {self.type}")

        # Initialise the tokens here and place on the device
        self.pad_token = (
            torch.tensor(self.number_tokeniser.encode(["<PAD>"]))
            .unsqueeze(0)
            .to(self.device)
        )
        self.sos_token = (
            torch.tensor(self.number_tokeniser.encode(["<SOS>"]))
            .unsqueeze(0)
            .to(self.device)
        )
        self.eos_token = (
            torch.tensor(self.number_tokeniser.encode(["<EOS>"]))
            .unsqueeze(0)
            .to(self.device)
        )
        self.skip_token = (
            torch.tensor(self.number_tokeniser.encode(["<SKIP>"]))
            .unsqueeze(0)
            .to(self.device)
        )
        self.x_token = (
            torch.tensor(self.number_tokeniser.encode(["X"]))
            .unsqueeze(0)
            .to(self.device)
        )

        self.model = self._load_model()

    def _load_model(self):
        model_loader = Loader(self.type, self.mode, self.device)
        return model_loader.model

    def __call__(self, tokenised_seqs: dict[str, TokenisedSequence], offsets):
        """
        This involves putting tokenised seqs into dataloader, making predictions,
        formating the output. Returning the numbered seqs to the user in the order in
        which they were given.
        """

        # NB: Provide a list of recommended batch sizes based on RAM and architecture

        dl = dataloader(self.batch_size, list(tokenised_seqs.values()))
        numbering = dict(zip(tokenised_seqs, self._predict_numbering(dl)))

        # Add offsets, where necessary.
        for key, value in offsets.items():
            # Catch long sequences which have an offset but fail numbering.
            if numbering[key]["query_start"] is None:
                continue
            numbering[key]["query_start"] += value
            numbering[key]["query_end"] += value

        return numbering

    def _predict_numbering(self, dl):
        """
        1 Runs the autoregressive inference loop which takes batches of sequences and
        predicts for the whole batch

        2 Converts token values to word and transfer to CPU  (if running on GPU)

        3 Works out the valid indicies - IMGT integer numbers predicted by the model

        4 find the EOS predicted by the model and the EOS of the input sequence

        5 iterates through each batch of predictions
            5A - check score is valid
            5B - if valid iterate over individual seq.
            5C - forward fill to end of the sequence for missed numbering
            5D - back fill for the start if the seq
            5E - Fill in up to 1 (starting IMGT residue) with gaps
            5F - Add gaps to nums where we are missing a number in the middle (27, 29)

        6 Populate the meta data dict and add to alignment

        Return the
        """
        if self.verbose:
            print(f"Making predictions on {len(dl)} batches.")

        numbering = []

        num = 0

        pad_token = self.pad_token
        sos_token = self.sos_token
        eos_token = self.eos_token
        skip_token = self.skip_token
        x_token = self.x_token

        ### 1 RUN AUTOREGRESSIVE INFERENCE LOOP OVER BATCHES

        with torch.no_grad():
            for X in dl:
                src = X.to(self.device)
                batch_size = src.shape[0]
                trg_len = src.shape[1] + 1  # Need to add 1 to include chain ID

                src_mask = self.model.make_src_mask(src)
                enc_src = self.model.encoder(src, src_mask)

                input = src[:, 0].unsqueeze(1)
                mask_input = src[:, 0].unsqueeze(1)
                cache = None

                max_input = torch.zeros(
                    batch_size, trg_len, device=self.device, dtype=torch.long
                )
                max_input[:, 0] = src[:, 0]

                scores = torch.zeros(
                    batch_size, trg_len - 1, device=self.device, dtype=torch.float
                )

                for t in range(1, trg_len):
                    trg_pad_mask, trg_causal_mask = self.model.make_trg_mask(mask_input)

                    output, cache = self.model.decoder(
                        input, enc_src, trg_pad_mask, trg_causal_mask, src_mask, cache
                    )

                    pred_token = output.argmax(2)[:, -1].unsqueeze(1)

                    scores[:, t - 1 : t] = output.topk(1, dim=2).values.squeeze(1)

                    max_input[:, t : t + 1] = pred_token

                    mask_input = max_input[:, : t + 1]

                    input = pred_token

                ### 2 tokenise and transfer the batch to cpu

                src_tokens = self.sequence_tokeniser.tokens[src.to("cpu")]
                pred_tokens = self.number_tokeniser.tokens[max_input.to("cpu")]

                ### 3 work out IMGT integer values predicted by model

                scores = scores.squeeze(
                    -1
                )  # Remove the last dim; shape becomes [batch_size, trg_len]

                mask = (
                    (max_input != skip_token)
                    & (max_input != x_token)
                    & (max_input != pad_token)
                    & (max_input != sos_token)
                )

                ### 4 Find the predicted end of sequence by model,
                # find actual end of input

                # Find first `True` (eos_token) along last dim (trg_len)
                eos_positions = max_input == eos_token
                # Get the indices (trg_len), for each batch
                first_eos_positions = torch.argmax(eos_positions.to(torch.int64), dim=1)

                # Same logic to find SRC EOS position
                src_eos_matrix = src == eos_token
                src_eos_positions = torch.argmax(src_eos_matrix.to(torch.int64), dim=1)

                # Check if no EOS token is found for each batch
                no_eos_found = ~(eos_positions.any(dim=1))
                # True if no EOS token is found in the row
                # Set the position to trg_len if no EOS is found

                first_eos_positions[no_eos_found] = torch.tensor(
                    trg_len - 1, device=self.device
                )

                ### 5   Iterate through each seq in batch

                for batch_no in range(batch_size):
                    error_occurred = False
                    num += 1
                    error_msg = None

                    # eos_position = first_eos_positions[batch_no]

                    # Code fix here. Ensure that if model has numbered beyond SRC EOS
                    # Then stop.
                    eos_position = min(
                        first_eos_positions[batch_no], src_eos_positions[batch_no] + 1
                    )

                    valid_indices = torch.arange(eos_position, device=self.device)[
                        mask[batch_no, :eos_position]
                    ]
                    valid_scores = scores[batch_no, valid_indices]

                    ### 5A   Check score is valid

                    if len(valid_indices) >= 50:
                        normalized_score = valid_scores.mean().item()
                    else:
                        normalized_score = 0.0
                        error_msg = "Less than 50 non insertion residues numbered."

                    if normalized_score < CUTOFF_SCORE:
                        numbering.append(
                            {
                                "numbering": None,
                                "chain_type": "F",
                                "score": normalized_score,
                                "query_start": None,
                                "query_end": None,
                                "error": error_msg or "Score less than cut off.",
                                "scheme": "imgt",
                            }
                        )
                        # skip the rest of the loop.
                        continue

                    ### 5B   Begin populating the numbering labels but iterating over
                    # each seq

                    residues, nums = [], []
                    backfill_residues = []

                    started = False
                    in_x_run, x_count = False, 0
                    start_index = None
                    end_index = None

                    # SRC is missing chain token + 1
                    src_eos_position = src_eos_positions[batch_no].item() + 1
                    eos_position = eos_position.item()

                    for seq_position in range(2, eos_position):
                        ###      Break at actual EOS  in the input sequence
                        if src_tokens[batch_no, seq_position - 1] == "<EOS>":
                            # The end index position in the sequence
                            # -3 is to accomodate the shifted register due to -
                            #  the <SOS>, chain token and python zero
                            end_index = seq_position - 3
                            break

                        ###      Break at SKIP tokens if numbering has started
                        elif (
                            pred_tokens[batch_no, seq_position] == "<SKIP>" and started
                        ):  # Break if hitting a skip post at the end.
                            end_index = seq_position - 3
                            break

                        ###      Work out when numbering begins, ignore SKIP tokens if
                        ###      not started
                        elif (
                            pred_tokens[batch_no, seq_position] == "<SKIP>"
                            and not started
                        ):  # Append as backfill up to the start.
                            backfill_residues.append(
                                src_tokens[batch_no, seq_position - 1]
                            )
                            continue

                        ###      If an instertion X is called, log as  in a run of X
                        elif pred_tokens[batch_no, seq_position] == "X":
                            x_count += 1
                            in_x_run = True

                        ###      If breaking out of a X run, construct the labels
                        elif (
                            isinstance(pred_tokens[batch_no, seq_position], int)
                            and in_x_run
                        ):
                            # This code breaks if we have a junk seq that
                            # has predicted runs of X (insertions)
                            # that are not bookended with integers
                            try:
                                construction = build_inward_list(
                                    length=x_count,
                                    # number before X began
                                    start_num=int(
                                        pred_tokens[
                                            batch_no, (seq_position - (x_count + 1))
                                        ]
                                    ),
                                    # current number
                                    end_num=int(pred_tokens[batch_no, seq_position]),
                                )

                                # Add the construction over the previous sequence
                                nums[(seq_position - x_count) : seq_position] = (
                                    construction
                                )
                                # add the end
                                nums.append(
                                    (int(pred_tokens[batch_no, seq_position]), " ")
                                )
                                in_x_run = False
                                x_count = 0

                            except ValueError as e:
                                # Capture the error message from the exception
                                captured_error = str(e)
                                numbering.append(
                                    {
                                        "numbering": None,
                                        "chain_type": "F",
                                        "score": normalized_score,
                                        "query_start": None,
                                        "query_end": None,
                                        "error": "Could not apply numbering: "
                                        f"{captured_error}",
                                        "scheme": "imgt",
                                    }
                                )
                                error_occurred = True
                                break

                        ###      No conditions have been found - it is a number label,
                        # append to nums
                        else:
                            try:
                                nums.append(
                                    (int(pred_tokens[batch_no, seq_position]), " ")
                                )
                            except ValueError as e:
                                # Capture the error message from the exception
                                captured_error = str(e)
                                numbering.append(
                                    {
                                        "numbering": None,
                                        "chain_type": "F",
                                        "score": normalized_score,
                                        "query_start": None,
                                        "query_end": None,
                                        "error": "Could not apply numbering: "
                                        f"{captured_error}",
                                        "scheme": "imgt",
                                    }
                                )
                                error_occurred = True
                                break

                        ###      After each iteration through the sequence append the
                        # sequence residue
                        residues.append(src_tokens[batch_no, seq_position - 1])

                        if not started:
                            start_index = seq_position - 2
                        started = True

                    if error_occurred:
                        continue

                    # Assign an end index before entering forwardfill
                    if not end_index:
                        end_index = eos_position - 3
                        # eos_position - 1: Moves to the token before <EOS>,
                        # excluding the <EOS> marker itself.
                        # Subtracting an additional 1 for SOS and 1 for CLS:
                        # Adjusts further to skip over these two tokens.

                    ## Check for duplicates
                    if len(nums) != len(set(nums)):
                        numbering.append(
                            {
                                "numbering": None,
                                "chain_type": "F",
                                "score": normalized_score,
                                "query_start": None,
                                "query_end": None,
                                "error": "Model predicted duplicate numbers",
                                "scheme": "imgt",
                            }
                        )
                        # break out of the loop
                        continue

                    ### 5C   Perform forward fill to end of the sequence, if
                    # missed numbering

                    ##  ANARCII sometimes doesn't continue numbering to end of seq
                    # Solution: Identify residues remaining after the EOS
                    # Decide forward fill to 127 (KL) /128 (H) needs to occur.

                    # The last number depends on chain type - check type here.
                    if pred_tokens[batch_no, 1] in ["H", "A", "G"]:
                        last_num = 128
                    else:
                        last_num = 127

                    try:
                        last_predicted_num = int(
                            pred_tokens[batch_no, eos_position - 1]
                        )
                    except ValueError:
                        last_predicted_num = last_num

                    ### DEBUG ONLY ###
                    # if src_tokens[batch_no, eos_position - 1] == "<PAD>":
                    #     print(src_tokens[batch_no, :])
                    #     print(pred_tokens[batch_no, :])

                    # print(src_tokens[batch_no, eos_position - 1])
                    # print(last_num, last_predicted_num)
                    # print(last_predicted_num != last_num)

                    if (
                        src_tokens[batch_no, eos_position - 1] not in ["<EOS>", "<PAD>"]
                        and last_predicted_num != last_num
                        and last_predicted_num > 119
                    ):
                        # How far is EOS from 128?
                        missing_count = last_num - int(
                            pred_tokens[batch_no, eos_position - 1]
                        )

                        # How much is left of the source to number?
                        seq_remainder = int(src_eos_position) - int(eos_position)

                        missing_end_nums = [
                            (x, " ")
                            for x in range(last_predicted_num + 1, last_num + 1)
                        ]

                        missing_end_nums = missing_end_nums[:seq_remainder]

                        # # DEBUG PURPOSE ONLY
                        # print("\n")
                        # for i in range(
                        #     eos_position-3,
                        #     eos_position + min(missing_count, seq_remainder)):
                        #     print(
                        #         "\t", src_tokens[batch_no, i + 0],
                        #         "\t", pred_tokens[batch_no, i + 1],
                        #     )

                        missing_end_residues = []
                        for i in range(
                            eos_position,
                            eos_position + min(missing_count, seq_remainder),
                        ):
                            missing_end_residues.append(src_tokens[batch_no, i - 1])

                        # print("Last:\t", last_num)
                        # print("Last pred num:\t", last_predicted_num)
                        # print("Missing count:\t", missing_count)

                        # print("Missing num:\t", missing_end_nums)
                        # print("Missing res:\t", missing_end_residues)

                        # # Append the misssing labels to seq and nums:
                        nums = nums + missing_end_nums
                        residues = residues + missing_end_residues

                        end_index = end_index + len(missing_end_nums)

                    ### 5D   Perform backfill for missed start of sequence, if missed
                    # numbering

                    # This step ensures that first and last nums will always be
                    # integers - does not proceed if not.
                    try:
                        first_num = int(nums[0][0])  # get first number
                        last_num = int(nums[-1][0])  # get last number
                    except (IndexError, ValueError) as e:
                        # When numbering has failed, `nums` is an empty list.
                        # For some non-antibody/TCR sequences that the model does not
                        # recognise, the first number can be a string, like an EOS or
                        # an X token. End the loop here and move on to the next seq
                        # in the batch.
                        captured_error = str(e)
                        numbering.append(
                            {
                                "numbering": None,
                                "chain_type": "F",
                                "score": normalized_score,
                                "query_start": None,
                                "query_end": None,
                                "error": f"Could not apply numbering: {captured_error}",
                                "scheme": "imgt",
                            }
                        )
                        continue

                    # Should not do this before 10 in case of failure to
                    # identify the gap.
                    if first_num > 1 and first_num < 9 and len(backfill_residues) > 0:
                        # This creates a list from 1 to first_num - 1
                        vals = list(range(1, first_num))
                        # the problem here is if there is a lot of junk...
                        vals = vals[-len(backfill_residues) :]
                        nums = [(i, " ") for i in vals] + nums
                        residues = list(backfill_residues[-len(vals) :]) + residues

                        # Adjust the start index for the backfill
                        start_index = start_index - len(vals)

                    ### 5E Fill in up to 1 (starting IMGT residue) with gaps
                    first_num = int(nums[0][0])  # get first number again - may change
                    for missing_num in range(
                        first_num - 1, 0, -1
                    ):  # Start from first_num - 1, stop at 1, step by -1
                        nums.insert(0, (missing_num, " "))
                        residues.insert(0, "-")

                    ### 5F Add gaps to nums where we are missing a number:
                    # e.g. predicted labels are 91 L, 93 K. convert to >>
                    # 91 L, 92 -, 93 K
                    i = 1
                    while i < len(nums):
                        if (int(nums[i][0]) - 1) > int(nums[i - 1][0]):
                            nums.insert(i, (int(nums[i - 1][0]) + 1, " "))
                            residues.insert(i, "-")
                        else:
                            i += 1  # Only increment if no insertion is made

                    # Ensure the last number is 128 >>>>>
                    last_num = int(nums[-1][0])
                    for missing_num in range(last_num + 1, 129):
                        nums.append((missing_num, " "))
                        residues.append("-")

                    ### 6 Populate the meta data dict and append to alignment list

                    # Successful - append.
                    numbering.append(
                        {
                            "numbering": list(zip(nums, residues)),
                            "chain_type": str(pred_tokens[batch_no, 1]),
                            "score": normalized_score,
                            "query_start": start_index,
                            "query_end": end_index,
                            "error": None,
                            "scheme": "imgt",
                        }
                    )

            return numbering
