import json
import os

import torch
import torch.nn.functional as F

from anarcii.classifii import model
from anarcii.inference.utils import dataloader
from anarcii.input_data_processing.tokeniser import Tokeniser

type_tokens = {"A": "antibody", "T": "tcr"}


class TypeTokeniser(Tokeniser):
    def __init__(self, vocab_type="protein"):
        self.vocab_type = vocab_type
        self.pad = "<PAD>"
        self.start = "<SOS>"
        self.end = "<EOS>"

        if self.vocab_type == "protein":
            self.vocab = [
                self.pad,
                self.start,
                self.end,
                *([x.upper() for x in "acdefghiklmnpqrstvwXy"]),
            ]

        elif self.vocab_type == "number":
            self.vocab = [self.pad, self.start, *type_tokens]
        else:
            raise ValueError(f"Vocab type {vocab_type} not supported")

        super().__init__()


class TypeLoader:
    def __init__(self, device):
        self.device = device
        self.script_dir = os.path.dirname(os.path.abspath(__file__))
        params = self._load_params()

        self.INPUT_DIM = params["INPUT_DIM"]
        self.OUTPUT_DIM = params["OUTPUT_DIM"]
        self.HID_DIM = params["HID_DIM"]
        self.ENC_LAYERS = self.DEC_LAYERS = params["LAYERS"]
        self.ENC_HEADS = self.DEC_HEADS = params["HEADS"]
        self.ENC_PF_DIM = self.DEC_PF_DIM = params["HID_DIM"] * 4
        self.ENC_DROPOUT = self.DEC_DROPOUT = params["DROPOUT"]
        self.SRC_PAD_IDX = params["SRC_PAD_IDX"]
        self.TRG_PAD_IDX = params["TRG_PAD_IDX"]

        self.model = self._load_model()

    def _load_params(self):
        params_path = os.path.join(self.script_dir, "params.json")

        with open(params_path) as file:
            params = json.load(file)
        return params

    def _load_model(self):
        ENC = model.Encoder(
            self.INPUT_DIM,
            self.HID_DIM,
            self.ENC_LAYERS,
            self.ENC_HEADS,
            self.ENC_PF_DIM,
            self.ENC_DROPOUT,
            self.device,
        )

        DEC = model.Decoder(
            self.OUTPUT_DIM,
            self.HID_DIM,
            self.DEC_LAYERS,
            self.DEC_HEADS,
            self.DEC_PF_DIM,
            self.DEC_DROPOUT,
            self.device,
        )

        S2S = model.S2S(ENC, DEC, self.SRC_PAD_IDX, self.TRG_PAD_IDX, self.device)

        model_path = os.path.join(self.script_dir, "classifii.pt")
        S2S.load_state_dict(
            torch.load(model_path, map_location=self.device, weights_only=True)
        )

        S2S.eval()
        return S2S


class Classifii:
    def __init__(self, batch_size, device):
        self.batch_size = batch_size
        self.device = device
        self.aa = TypeTokeniser("protein")
        self.num = TypeTokeniser("number")
        self.model = TypeLoader(self.device).model

    def __call__(self, sequences: dict[str, str]) -> dict[str, dict[str, str]]:
        tokenized_seqs = []
        # Capped at 235 for now.
        for seq in sequences.values():
            bookend_seq = [self.aa.start, *seq[:235], self.aa.end]
            try:
                tokenized_seqs.append(torch.from_numpy(self.aa.encode(bookend_seq)))
            except KeyError as e:
                print(
                    f"Sequence could not be numbered. Contains an invalid residue: {e}"
                )
                tokenized_seqs.append(torch.from_numpy(self.aa.encode(["F"])))

        dl = dataloader(self.batch_size, tokenized_seqs)
        classes = self._classify(dl)

        grouped_sequences = {type_tokens[key]: {} for key in set(classes)}
        for classification, (name, sequence) in zip(classes, sequences.items()):
            grouped_sequences[type_tokens[classification]][name] = sequence

        return grouped_sequences

    def _classify(self, dl):
        preds = []
        with torch.no_grad():
            for X in dl:
                src = X.to(self.device)
                src_mask = self.model.make_src_mask(src)
                enc_src = self.model.encoder(src, src_mask)

                input = src[:, 0].unsqueeze(1)
                trg_pad_mask, trg_causal_mask = self.model.make_trg_mask(input)
                output = self.model.decoder(
                    input, enc_src, trg_pad_mask, trg_causal_mask, src_mask
                )

                probs = F.softmax(output, dim=-1)
                pred_token = probs.argmax(2)[:, -1].unsqueeze(1)
                preds += [x[0] for x in self.num.tokens[pred_token.to("cpu")]]
        return preds
