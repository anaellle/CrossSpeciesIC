import numpy as np
import torch

from distances import CosineSimilarity
from utils import common_functions as c_f
from utils import loss_and_miner_utils as lmu
from losses.generic_pair_loss import GenericPairLoss


class TupletMarginLoss(GenericPairLoss):
    def __init__(self, margin=5.73, scale=64, **kwargs):
        super().__init__(mat_based_loss=False, **kwargs)
        c_f.assert_distance_type(self, CosineSimilarity)
        self.margin = np.radians(margin)
        self.scale = scale
        self.add_to_recordable_attributes(
            list_of_names=["margin", "scale"], is_stat=False
        )
        self.add_to_recordable_attributes(
            list_of_names=["avg_pos_angle", "avg_neg_angle"], is_stat=True
        )

    # pos_pairs and neg_pairs already represent cos(theta)
    def _compute_loss(self, pos_pairs, neg_pairs, indices_tuple):
        a1, p, a2, _ = indices_tuple

        if len(a1) > 0 and len(a2) > 0:
            pos_angles = torch.acos(pos_pairs)
            self.set_stats(pos_angles, neg_pairs)
            pos_pairs = torch.cos(pos_angles - self.margin)
            pos_pairs = pos_pairs.unsqueeze(1)
            neg_pairs = neg_pairs.repeat(pos_pairs.size(0), 1)
            inside_exp = self.scale * (neg_pairs - pos_pairs)
            keep_mask = a2.unsqueeze(0) == a1.unsqueeze(1)
            loss = lmu.logsumexp(inside_exp, keep_mask=keep_mask, add_one=True, dim=1)
            return {
                "loss": {
                    "losses": loss,
                    "indices": (a1, p),
                    "reduction_type": "pos_pair",
                }
            }
        return self.zero_losses()

    def get_default_distance(self):
        return CosineSimilarity()

    def set_stats(self, pos_angles, neg_pairs):
        if self.collect_stats:
            with torch.no_grad():
                neg_angles = torch.acos(neg_pairs)
                self.avg_pos_angle = np.degrees(torch.mean(pos_angles).item())
                self.avg_neg_angle = np.degrees(torch.mean(neg_angles).item())
