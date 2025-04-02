import math
import torch
import torch.nn as nn

class PositionalEncoding(nn.Module):
    def __init__(self, dim_model, dropout=0.1, max_len=5000):
        super(PositionalEncoding, self).__init__()
        self.dropout = nn.Dropout(p=dropout)
        position = torch.arange(0, max_len, dtype=torch.float32).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, dim_model, 2).float() * (-math.log(10000.0) / dim_model)
        )
        pe = torch.zeros(1, max_len, dim_model)
        pe[0, :, 0::2] = torch.sin(position * div_term)
        pe[0, :, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe)

    def forward(self, x):
        x = x + self.pe[:, : x.size(1), :]
        return self.dropout(x)


class TransformerModel(nn.Module):
    def __init__(
        self,
        input_dim,
        num_heads=8,
        num_layers=1,
        dim_feedforward=128,
        dropout=0.0
    ):
        super(TransformerModel, self).__init__()
        self.input_fc = nn.Linear(input_dim, dim_feedforward)
        self.pos_encoder = PositionalEncoding(dim_feedforward, dropout)
        encoder_layer = nn.TransformerEncoderLayer(
            d_model=dim_feedforward,
            nhead=num_heads,
            dim_feedforward=dim_feedforward * 4,
            dropout=dropout,
            activation='relu',
            batch_first=True
        )
        self.transformer_encoder = nn.TransformerEncoder(encoder_layer, num_layers=num_layers)

    def forward(self, src):
        src = self.input_fc(src)
        if src.dim() == 2:
            src = src.unsqueeze(1)
        src = self.pos_encoder(src)
        output = self.transformer_encoder(src)
        output = output.mean(dim=1)
        return output


class CombinedModel(nn.Module):
    def __init__(self, model1_params, combined_dim, output_dim):
        super(CombinedModel, self).__init__()
        self.model1 = TransformerModel(
            input_dim=model1_params['input_dim'],
            num_heads=model1_params['num_heads'],
            num_layers=model1_params['num_layers'],
            dim_feedforward=model1_params['dim_feedforward'],
            dropout=model1_params['dropout']
        )
        feature_dim = model1_params['dim_feedforward']
        self.fc = nn.Linear(feature_dim, combined_dim)
        self.relu = nn.ReLU()
        self.output_layer = nn.Linear(combined_dim, output_dim)

    def forward(self, input1):
        output1 = self.model1(input1)
        x = self.relu(self.fc(output1))
        return self.output_layer(x)
