import torch
import torch.nn as nn

seq_max_len = 310


class EncoderLayer(nn.Module):
    def __init__(self, hid_dim, n_heads, pf_dim, dropout, device):
        super().__init__()

        self.self_attn_layer_norm = nn.LayerNorm(hid_dim)
        self.ff_layer_norm = nn.LayerNorm(hid_dim)
        self.self_attention = EncoderMultiHeadAttentionLayer(
            hid_dim, n_heads, dropout, device
        )
        self.positionwise_feedforward = PositionwiseFeedforwardLayer(
            hid_dim, pf_dim, dropout
        )
        self.dropout = nn.Dropout(dropout)

    def forward(self, src, src_mask):
        # src = [batch size, src len, hid dim]
        # src_mask = [batch size, src len]

        # self attention
        _src, _ = self.self_attention(src, src, src, src_mask)

        # dropout, residual connection and layer norm
        src = self.self_attn_layer_norm(src + self.dropout(_src))
        # src = [batch size, src len, hid dim]

        # positionwise feedforward
        _src = self.positionwise_feedforward(src)

        # dropout, residual and layer norm
        src = self.ff_layer_norm(src + self.dropout(_src))
        # src = [batch size, src len, hid dim]

        return src


class Encoder(nn.Module):
    def __init__(
        self,
        input_dim,
        hid_dim,
        n_layers,
        n_heads,
        pf_dim,
        dropout,
        device,
        max_length=seq_max_len,
    ):
        super().__init__()

        self.device = device
        self.tok_embedding = nn.Embedding(input_dim, hid_dim)
        self.pos_embedding = nn.Embedding(max_length, hid_dim)
        self.layers = nn.ModuleList(
            [
                EncoderLayer(hid_dim, n_heads, pf_dim, dropout, device)
                for _ in range(n_layers)
            ]
        )
        self.dropout = nn.Dropout(dropout)
        self.scale = torch.sqrt(torch.FloatTensor([hid_dim])).to(device)

    def forward(self, src, src_mask):
        # src = [batch size, src len]
        # src_mask = [batch size, src len]
        batch_size = src.shape[0]
        src_len = src.shape[1]

        pos = (
            torch.arange(0, src_len).unsqueeze(0).repeat(batch_size, 1).to(self.device)
        )
        # pos = [batch size, src len]

        src = self.dropout(
            (self.tok_embedding(src) * self.scale) + self.pos_embedding(pos)
        )
        # src = [batch size, src len, hid dim]

        for layer in self.layers:
            src = layer(src, src_mask)

        # src = [batch size, src len, hid dim]
        return src


class EncoderMultiHeadAttentionLayer(nn.Module):
    def __init__(self, hid_dim, n_heads, dropout, device):
        super().__init__()

        # Initialize the built-in multi-head attention layer
        self.multihead_attn = nn.MultiheadAttention(
            hid_dim, n_heads, dropout=dropout, batch_first=True
        )

        # Ensure the multi-head attention layer and additional layers are moved to the
        # correct device
        self.device = device
        self.to(device)

    def forward(self, query, key, value, mask):
        # Forward pass through the built-in MultiheadAttention layer
        attn_output, attn_output_weights = self.multihead_attn(
            query, key, value, key_padding_mask=mask
        )

        return attn_output, attn_output_weights


class DecoderMultiHeadAttentionLayer(nn.Module):
    def __init__(self, hid_dim, n_heads, dropout, device):
        super().__init__()

        # Initialize the built-in multi-head attention layer
        self.multihead_attn = nn.MultiheadAttention(
            hid_dim, n_heads, dropout=dropout, batch_first=True
        )

        # Ensure the multi-head attention layer and additional layers are moved to the
        # correct device
        self.device = device
        self.to(device)

    def forward(self, query, key, value, trg_pad_mask, trg_causal_mask):
        # Forward pass through the built-in MultiheadAttention layer
        attn_output, attn_output_weights = self.multihead_attn(
            query,
            key,
            value,
            key_padding_mask=trg_pad_mask,
            attn_mask=trg_causal_mask,
            is_causal=True,
        )
        return attn_output, attn_output_weights


class PositionwiseFeedforwardLayer(nn.Module):
    def __init__(self, hid_dim, pf_dim, dropout):
        super().__init__()

        self.fc_1 = nn.Linear(hid_dim, pf_dim)
        self.fc_2 = nn.Linear(pf_dim, hid_dim)

        self.dropout = nn.Dropout(dropout)

    def forward(self, x):
        # x = [batch size, seq len, hid dim]

        x = self.dropout(torch.relu(self.fc_1(x)))
        # x = [batch size, seq len, pf dim]

        x = self.fc_2(x)
        # x = [batch size, seq len, hid dim]

        return x


class DecoderLayer(nn.Module):
    def __init__(self, hid_dim, n_heads, pf_dim, dropout, device):
        super().__init__()

        self.self_attn_layer_norm = nn.LayerNorm(hid_dim)
        self.enc_attn_layer_norm = nn.LayerNorm(hid_dim)
        self.ff_layer_norm = nn.LayerNorm(hid_dim)
        self.self_attention = DecoderMultiHeadAttentionLayer(
            hid_dim, n_heads, dropout, device
        )
        self.encoder_attention = EncoderMultiHeadAttentionLayer(
            hid_dim, n_heads, dropout, device
        )
        self.positionwise_feedforward = PositionwiseFeedforwardLayer(
            hid_dim, pf_dim, dropout
        )
        self.dropout = nn.Dropout(dropout)

    def forward(self, trg, enc_src, trg_pad_mask, trg_causal_mask, src_mask):
        # trg = [batch size, trg len, hid dim]
        # enc_src = [batch size, src len, hid dim]

        # src_mask = [batch size, 1, 1, src len]

        # self attention
        _trg, _ = self.self_attention(trg, trg, trg, trg_pad_mask, trg_causal_mask)

        # dropout, residual connection and layer norm
        trg = self.self_attn_layer_norm(trg + self.dropout(_trg))
        # trg = [batch size, trg len, hid dim]

        # encoder attention
        _trg, attention = self.encoder_attention(trg, enc_src, enc_src, src_mask)

        # dropout, residual connection and layer norm
        trg = self.enc_attn_layer_norm(trg + self.dropout(_trg))
        # trg = [batch size, trg len, hid dim]

        # positionwise feedforward
        _trg = self.positionwise_feedforward(trg)

        # dropout, residual and layer norm
        trg = self.ff_layer_norm(trg + self.dropout(_trg))

        # trg = [batch size, trg len, hid dim]
        # attention = [batch size, n heads, trg len, src len]

        return trg, attention


class Decoder(nn.Module):
    def __init__(
        self,
        output_dim,
        hid_dim,
        n_layers,
        n_heads,
        pf_dim,
        dropout,
        device,
        max_length=seq_max_len,
    ):
        super().__init__()

        self.device = device
        self.output_dim = output_dim
        self.tok_embedding = nn.Embedding(output_dim, hid_dim)
        self.pos_embedding = nn.Embedding(max_length, hid_dim)

        self.layers = nn.ModuleList(
            [
                DecoderLayer(hid_dim, n_heads, pf_dim, dropout, device)
                for _ in range(n_layers)
            ]
        )

        self.fc_out = nn.Linear(hid_dim, output_dim)
        self.dropout = nn.Dropout(dropout)
        self.scale = torch.sqrt(torch.FloatTensor([hid_dim])).to(device)

    def forward(self, trg, enc_src, trg_pad_mask, trg_causal_mask, src_mask):
        # trg = [batch size, trg len]
        # enc_src = [batch size, src len, hid dim]
        # src_mask = [batch size, src len]

        # print(
        #     "Dec shapes: ",
        #     trg.shape,
        #     enc_src.shape,
        #     trg_pad_mask.shape,
        #     trg_causal_mask.shape,
        #     src_mask.shape,
        # )

        batch_size = trg.shape[0]
        trg_len = trg.shape[1]
        pos = (
            torch.arange(0, trg_len, device=self.device)
            .unsqueeze(0)
            .repeat(batch_size, 1)
        )
        # pos = [batch size, trg len]

        trg = self.dropout(
            (self.tok_embedding(trg) * self.scale) + self.pos_embedding(pos)
        )
        # trg = [batch size, trg len, hid dim]

        for layer in self.layers:
            trg, _ = layer(trg, enc_src, trg_pad_mask, trg_causal_mask, src_mask)

        # trg = [batch size, trg len, hid dim]
        # attention = [batch size, n heads, trg len, src len]

        output = self.fc_out(trg)
        # output = [batch size, trg len, output dim]

        return output


class S2S(nn.Module):
    def __init__(self, encoder, decoder, src_pad_idx, trg_pad_idx, device):
        super().__init__()

        self.encoder = encoder
        self.decoder = decoder
        self.src_pad_idx = src_pad_idx
        self.trg_pad_idx = trg_pad_idx
        self.device = device
        self.to(device)

    def make_src_mask(self, src):
        # src = [batch size, src len]

        src_mask = (src == self.src_pad_idx).to(self.device)
        # src_mask = [batch size, src len]

        return src_mask

    def make_trg_mask(self, trg):
        # trg = [batch size, trg len]

        # Create a padding mask of shape [batch size, 1, trg len]
        trg_pad_mask = (trg == self.trg_pad_idx).to(self.device)

        # trg_pad_mask = [batch size, trg len]
        # Adjusting shape for broadcasting by adding an extra dimension for trg_len to
        # match causal_mask
        trg_len = trg.shape[1]

        # Create a subsequence mask of shape [trg len, trg len]
        trg_causal_mask = ~torch.tril(
            torch.ones((trg_len, trg_len), device=self.device)
        ).bool()
        # trg_causal_mask = [trg len, trg len]

        return trg_pad_mask, trg_causal_mask

    def forward(self, src, trg):
        src_mask = self.make_src_mask(src)
        trg_pad_mask, trg_causal_mask = self.make_trg_mask(trg)

        enc_src = self.encoder(src, src_mask)
        output = self.decoder(trg, enc_src, trg_pad_mask, trg_causal_mask, src_mask)

        return output
