import torch

from anarcii.inference.window_selector import WindowFinder


def split_seq(seq: str, n_jump: int, window_size: int = 90) -> list[str]:
    jump = n_jump
    num = (len(seq) - window_size) // jump
    ls = [seq[(jump * x) : (jump * x + window_size)] for x in range(num)]
    return ls


def pick_windows(
    seqs: list[str], model: WindowFinder, fallback: bool = False
) -> list[int] | int:
    # Find the index of the highest scoring window
    aa = model.sequence_tokeniser
    ls = []

    for seq in seqs:
        bookend_seq = [aa.start] + list(seq) + [aa.end]
        try:
            tokenised_seq = torch.from_numpy(aa.encode(bookend_seq))
            ls.append(tokenised_seq)
        except KeyError as e:
            print(f"Sequence could not be numbered. Contains an invalid residue: {e}")
            ls.append([])

    max_index = model(ls, fallback)
    return max_index
