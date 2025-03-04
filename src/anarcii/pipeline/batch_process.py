from anarcii.input_data_processing.sequences import SequenceProcessor


def batch_process(ls, model, window_model, verbose, txt_file):
    counter = 1
    for chunk in ls:
        print(f"\nChunk: {counter} of {len(ls)}.")

        sequences = SequenceProcessor(chunk, model, window_model, verbose)
        processed_seqs, offsets = sequences.process_sequences()

        # Process and write to the temp file
        numbered_seqs = model(processed_seqs, offsets)

        with open(txt_file, "a") as file:
            for item in numbered_seqs:
                file.write(repr(item) + "\n")

        counter += 1

    return []
