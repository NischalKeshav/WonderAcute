import numpy as np
import pytest

from qarch import qram
from qarch.passes import compress


@pytest.mark.parametrize("k", [1, 2, 3, 4, 5])
@pytest.mark.parametrize("d", [1, 3])
def test_table_lookup_reads_every_address(k, d):
    N = 2 ** k
    data = [int(v) for v in np.random.default_rng(k * 10 + d).integers(0, 2 ** d, N)]
    lk = qram.table_lookup(data, d)
    assert qram.check_lookup(lk, data)
    c = qram.table_lookup_counts(N, d)
    assert [i.operation.name for i in lk.circuit.data].count("ccx") == c["toffoli"]
    names = [i.operation.name for i in qram.table_lookup(data, d, exact=True).circuit.data]
    # 7 T gates per Toffoli in the exact decomposition
    assert sum(n in ("t", "tdg") for n in names) == c["t_count"]
    assert lk.circuit.num_qubits == c["qubits"]


@pytest.mark.parametrize("k", [1, 2, 3, 4])
@pytest.mark.parametrize("d", [1, 2])
def test_bucket_brigade_reads_every_address(k, d):
    N = 2 ** k
    data = [int(v) for v in np.random.default_rng(k * 7 + d).integers(0, 2 ** d, N)]
    lk = qram.bucket_brigade(data, d)
    assert qram.check_lookup(lk, data)
    c = qram.bucket_brigade_counts(N, d)
    names = [i.operation.name for i in lk.circuit.data]
    assert names.count("cswap") == c["cswap"]
    assert names.count("x") - sum(bin(v).count("1") for v in data) == c["x_absorbable"]
    assert lk.circuit.num_qubits == c["qubits"]


def test_pass_absorbs_router_x_gates():
    data = [1, 0, 1, 1, 0, 0, 1, 0]
    lk = qram.bucket_brigade(data)
    res = compress(lk.circuit)
    assert res.removed["x"] >= qram.bucket_brigade_counts(8)["x_absorbable"]
    # CSWAP is not diagonal, so the X on the router has to be flushed before it:
    # the honest saving is what remains after flushes
    assert res.flushed <= res.removed["x"]
