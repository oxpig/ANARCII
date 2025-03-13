from anarcii import Anarcii

model = Anarcii(
    seq_type="unknown",
    batch_size=128,
    cpu=False,
    ncpu=4,
    mode="accuracy",
    verbose=True,
    # max_seqs_len=5, # Ensure we test the write to output functionality.
)
model.number("../data/raw_data/unknown.fa")

model.to_text("../data/expected_data/unknown_expected_1.txt")
model.to_json("../data/expected_data/unknown_expected_1.json")
