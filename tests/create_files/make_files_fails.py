from anarcii import Anarcii

model = Anarcii(
    seq_type="unknown",
    batch_size=32,
    cpu=False,
    ncpu=12,
    mode="speed",
    verbose=True,
)
model.number("../data/raw_data/fails.fa")

model.to_scheme("imgt")
model.to_msgpack("../data/expected_data/fails_imgt_expected.msgpack")

model.to_scheme("kabat")
model.to_msgpack("../data/expected_data/fails_kabat_expected.msgpack")
