# tests/test_mirpy_bamnostic.py
from mirpy import tester

def test_mirpytest_runs(monkeypatch):
    called = {}
    def fake_alignmentfile(path, mode):
        called['path'] = path
        class Dummy:
            def fetch(self, **kw): return []
            def close(self): pass
        return Dummy()
    monkeypatch.setattr(tester.bn, "AlignmentFile", fake_alignmentfile)
    rc = tester.miRpyTest(["fake.bam"], n=3)
    assert rc == 0
    assert called['path'] == "fake.bam"
