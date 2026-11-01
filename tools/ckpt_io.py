#!/usr/bin/env python3
"""ckpt_io.py -- read and WRITE ngspice_augmentation's `savestate`/`loadstate`
checkpoint file format (card #512, envelope-following round 10), in pure
Python, with zero new ngspice C code.

`docs/checkpoint.md` (card #472) already characterises what the format
holds and confirms `loadstate` restores it byte-for-byte. This module goes
one step further: it parses the file into named fields and can re-serialise
a MODIFIED copy that `loadstate` will accept -- the missing piece round 9's
own `tools/tran_from_pss.py` flagged ("loadstate restores the saved state
exactly; you would need a way to write a modified state") and this card's
own brief asks to resolve by "check[ing] the savestate file format".

Format, read directly from `src/frontend/com_checkpoint.c`
(`com_savestate`/`com_loadstate`, ngspice-src-env's `aug` branch, unchanged
since patches 0032-0033 -- see docs/checkpoint.md): a sequence of RAW
`fwrite`/`fread` calls, no struct padding between them (each field is its
own read/write of `sizeof(*(ptr)) * n` bytes), native type sizes. On the
build this project runs (x86-64 Linux, gcc), the C types involved are:

    magic            8 bytes literal b"NGCKPT01"
    maxEqNum         int32          (header: topology signature)
    matSize          int32          (header)
    numStates        int32          (header)
    maxOrder         int32          (header)
    breakSize        int32          (header)
    CKTorder         int32          (scalar integration state)
    CKTmode          int64          (cktdefs.h: `long CKTmode` -- NOT int;
                                      found by reading cktdefs.h directly,
                                      not assumed -- a naive int32 guess here
                                      would silently misalign every field
                                      after it)
    CKTtime          float64
    CKTdelta         float64
    CKTsaveDelta     float64
    CKTdelmin        float64
    CKTminBreak      float64
    CKTdeltaOld[7]   7 x float64
    CKTbreaks[]      breakSize x float64      (only present if breakSize>0)
    rhsOld_present   int32
    CKTrhsOld[]      (matSize+1) x float64    (only if rhsOld_present)
    irhsOld_present  int32
    CKTirhsOld[]     (matSize+1) x float64    (only if irhsOld_present)
    8x:
      state_present  int32
      CKTstates[i][]  numStates x float64     (only if state_present)

Verified (not assumed) by `tests/test_envelope.py`'s
`test_roundtrip_byte_identical()`: read a real checkpoint written by
`savestate`, write it back unmodified, and diff the two files byte-for-byte
-- 0 bytes different, on every fixture this module is used with.

This module intentionally does NOT try to interpret `CKTstates[]` (the
per-device integration history -- charge/flux/current history in a layout
that is private to each device model's own `.._STATE_INDEX` constants, not
documented anywhere at the file-format level). `envelope_tran.py` reads it
back and rewrites it VERBATIM (never modified), and instead resets the
integration order/breakpoint table so the engine re-derives a consistent
history over the next few real steps -- see that module's own docstring
for why (the same "known state discontinuity -> low order restart"
mechanism `docs/checkpoint.md`'s own patches 0200/0201 already rely on for
the ordinary breakpoint case).
"""
import struct

MAGIC = b"NGCKPT01"

_I = "=i"   # C int, 4 bytes, standard size, native byte order
_Q = "=q"   # C long on x86-64 Linux, 8 bytes
_D = "=d"   # C double, 8 bytes


class CheckpointFormatError(RuntimeError):
    pass


def read_checkpoint(path):
    """Parse a `savestate`-written file into a dict of named fields.
    Raises CheckpointFormatError if the magic does not match (same check
    `com_loadstate` itself makes)."""
    with open(path, "rb") as fh:
        data = fh.read()
    off = 0

    def take(fmt, n=1):
        nonlocal off
        elt_size = struct.calcsize(fmt)
        size = elt_size * n
        if off + size > len(data):
            raise CheckpointFormatError(
                f"{path}: truncated at offset {off} reading {n} x {fmt!r} "
                f"(file is {len(data)} bytes)")
        vals = struct.unpack_from("=" + fmt[-1] * n, data, off)
        off += size
        return list(vals)

    magic = data[0:8]
    off = 8
    if magic != MAGIC:
        raise CheckpointFormatError(f"{path}: bad magic {magic!r}, not an ngspice checkpoint")

    maxEqNum, matSize, numStates, maxOrder, breakSize = take(_I, 5)
    ckt_order, = take(_I, 1)
    ckt_mode, = take(_Q, 1)
    (ckt_time, ckt_delta, ckt_saveDelta, ckt_delmin, ckt_minBreak) = take(_D, 5)
    deltaOld = take(_D, 7)

    breaks = take(_D, breakSize) if breakSize > 0 else []

    rhs_present, = take(_I, 1)
    rhsOld = take(_D, matSize + 1) if rhs_present else None

    irhs_present, = take(_I, 1)
    irhsOld = take(_D, matSize + 1) if irhs_present else None

    states = []
    for _ in range(8):
        present, = take(_I, 1)
        if present and numStates > 0:
            states.append(take(_D, numStates))
        elif present:
            states.append([])
        else:
            states.append(None)

    if off != len(data):
        raise CheckpointFormatError(
            f"{path}: {len(data) - off} trailing bytes not accounted for "
            "(format assumption is wrong for this build)")

    return dict(
        maxEqNum=maxEqNum, matSize=matSize, numStates=numStates,
        maxOrder=maxOrder, breakSize=breakSize,
        CKTorder=ckt_order, CKTmode=ckt_mode, CKTtime=ckt_time,
        CKTdelta=ckt_delta, CKTsaveDelta=ckt_saveDelta, CKTdelmin=ckt_delmin,
        CKTminBreak=ckt_minBreak, CKTdeltaOld=deltaOld, CKTbreaks=breaks,
        CKTrhsOld=rhsOld, CKTirhsOld=irhsOld, CKTstates=states,
    )


def write_checkpoint(path, st):
    """Serialise a dict in the same shape `read_checkpoint` returns back to
    the on-disk format `loadstate` reads. Field-for-field inverse of
    `read_checkpoint` -- an unmodified read/write round trip is
    byte-identical to the original file (tests/test_envelope.py checks
    this directly)."""
    buf = bytearray()
    buf += MAGIC
    buf += struct.pack("=5i", st["maxEqNum"], st["matSize"], st["numStates"],
                        st["maxOrder"], st["breakSize"])
    buf += struct.pack(_I, st["CKTorder"])
    buf += struct.pack(_Q, st["CKTmode"])
    buf += struct.pack("=5d", st["CKTtime"], st["CKTdelta"], st["CKTsaveDelta"],
                        st["CKTdelmin"], st["CKTminBreak"])
    deltaOld = list(st["CKTdeltaOld"])
    if len(deltaOld) != 7:
        raise CheckpointFormatError(f"CKTdeltaOld must have 7 entries, got {len(deltaOld)}")
    buf += struct.pack("=7d", *deltaOld)

    breaks = st["CKTbreaks"]
    if st["breakSize"] > 0:
        if len(breaks) != st["breakSize"]:
            raise CheckpointFormatError(
                f"breakSize={st['breakSize']} but CKTbreaks has {len(breaks)} entries")
        buf += struct.pack(f"={st['breakSize']}d", *breaks)

    rhsOld = st["CKTrhsOld"]
    buf += struct.pack(_I, 1 if rhsOld is not None else 0)
    if rhsOld is not None:
        n = st["matSize"] + 1
        if len(rhsOld) != n:
            raise CheckpointFormatError(f"CKTrhsOld must have {n} entries, got {len(rhsOld)}")
        buf += struct.pack(f"={n}d", *rhsOld)

    irhsOld = st["CKTirhsOld"]
    buf += struct.pack(_I, 1 if irhsOld is not None else 0)
    if irhsOld is not None:
        n = st["matSize"] + 1
        if len(irhsOld) != n:
            raise CheckpointFormatError(f"CKTirhsOld must have {n} entries, got {len(irhsOld)}")
        buf += struct.pack(f"={n}d", *irhsOld)

    states = st["CKTstates"]
    if len(states) != 8:
        raise CheckpointFormatError(f"CKTstates must have 8 entries, got {len(states)}")
    for arr in states:
        if arr is None:
            buf += struct.pack(_I, 0)
        else:
            buf += struct.pack(_I, 1)
            n = st["numStates"]
            if n > 0:
                if len(arr) != n:
                    raise CheckpointFormatError(f"a CKTstates[] entry must have {n} doubles, got {len(arr)}")
                buf += struct.pack(f"={n}d", *arr)

    with open(path, "wb") as fh:
        fh.write(bytes(buf))


if __name__ == "__main__":
    import sys
    st = read_checkpoint(sys.argv[1])
    print(f"maxEqNum={st['maxEqNum']} matSize={st['matSize']} numStates={st['numStates']} "
          f"maxOrder={st['maxOrder']} breakSize={st['breakSize']}")
    print(f"CKTorder={st['CKTorder']} CKTmode={st['CKTmode']} CKTtime={st['CKTtime']!r} "
          f"CKTdelta={st['CKTdelta']!r} CKTsaveDelta={st['CKTsaveDelta']!r}")
    print(f"CKTbreaks[:5]={st['CKTbreaks'][:5]!r} ({len(st['CKTbreaks'])} total)")
    print(f"CKTrhsOld len={len(st['CKTrhsOld']) if st['CKTrhsOld'] else None}")
    print(f"CKTstates present: {[a is not None for a in st['CKTstates']]}")
