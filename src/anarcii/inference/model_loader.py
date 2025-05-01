import importlib.resources as pkg_resources
import json
import os

import torch

from . import model


class Loader:
    def __init__(self, sequence_type, mode, device):
        self.device = device
        self.type = sequence_type
        self.mode = mode

        # Based on the user inputs this loads the model parameters
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

        self.model_path = self._get_model_path()
        self.model = self._load_model()

    def _load_params(self):
        if self.type == "shark":
            param_filename = f"{self.type}_4_2_128_512.json"
        elif self.mode == "speed":
            param_filename = f"{self.type}_4_1_128_512.json"
        elif self.mode == "accuracy":
            param_filename = f"{self.type}_4_2_128_512.json"
        else:
            raise ValueError(
                "Invalid mode specified. Choose either 'speed' or 'accuracy' or "
                "'shark' (aliases 'vnar', 'vhh')."
            )

        param_path = pkg_resources.files("anarcii.models").joinpath(
            self.type, param_filename
        )

        if not param_path.exists():
            raise FileNotFoundError(f"Parameter file not found: {param_path}")

        with param_path.open("r") as file:
            params = json.load(file)

        return params

    def _get_model_path(self):
        model_filename = f"{self.type}_{self.ENC_HEADS}_{self.ENC_LAYERS}_{self.HID_DIM}_{self.ENC_PF_DIM}.pt"  # noqa: E501

        model_path = pkg_resources.files("anarcii.models").joinpath(
            self.type, model_filename
        )

        if not model_path.exists():
            raise FileNotFoundError(f"Model file not found: {model_path}")

        return str(model_path)

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

        if not os.path.exists(self.model_path):
            raise FileNotFoundError(f"Model file not found: {self.model_path}")

        S2S.load_state_dict(
            torch.load(self.model_path, map_location=self.device, weights_only=True)
        )

        S2S.eval()

        return S2S
