from anarcii import Anarcii

model = Anarcii(seq_type="antibody", mode="speed")
model.number("../data/raw_data/unknown.fa")
model.to_json("../data/expected_data/cli_a.json")


model = Anarcii(seq_type="tcr", mode="speed")
model.number("../data/raw_data/unknown.fa")
model.to_json("../data/expected_data/cli_b.json")


model = Anarcii(seq_type="shark", mode="speed")
model.number("../data/raw_data/unknown.fa")
model.to_json("../data/expected_data/cli_c.json")


model = Anarcii(seq_type="unknown", mode="speed")
model.number("../data/raw_data/unknown.fa")
model.to_json("../data/expected_data/cli_d.json")
