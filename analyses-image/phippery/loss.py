import torch


def mse(y_true, y_predicted):
    return torch.nn.functional.mse_loss(y_true.squeeze(), y_predicted.squeeze())
