"""Replay the Windows Gmsh 4.15.2 power oracle (no Tessella imports).

Requires Python gmsh==4.15.2 on Windows x86_64. Run normally to verify the
committed independent corpus. --regenerate explicitly replaces corpus/metadata.
--gmsh-executable PATH also checks the standalone parser at PC64, nearest mode.
"""
from pathlib import Path
import argparse
import collections
import ctypes
import hashlib
import json
import math
import os
import platform
import random
import struct
import subprocess
import tempfile
import gmsh

HERE=Path(__file__).resolve().parent
MODES=[(precision|rounding) for precision in (0x0200,0x0300,0x0000)
       for rounding in (0,0x0400,0x0800,0x0c00)]

def bits(value):
    return int.from_bytes(struct.pack("<d",value),"little")

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def case_catalog():
    cases = []
    def add(group, x, y):
        cases.append((group, float(x), float(y)))
    for x in (math.nextafter(1., 2.), math.nextafter(1., 0.)):
        for negative in (False, True):
            for n in (-2**31-1, -2**31, -2**31+1, 2**31-1, 2**31,
                      2**31+1, -2**53+1, 2**53-1, 2**53, 2**53+2, 2**62):
                add("signed_int32_cutoff", -x if negative else x, n)
    rng = random.Random(221097)
    for i in range(2000):
        if i % 4 == 0:
            x = math.ldexp(1 + rng.random(), rng.randrange(-1074, 1024))
        elif i % 4 == 1:
            x = 0.5 + rng.random()
        elif i % 4 == 2:
            x = math.nextafter(1., 2. if rng.randrange(2) else 0.)
        else:
            x = math.ldexp(1 + rng.random(), rng.randrange(510, 540))
        x *= -1 if rng.randrange(2) else 1
        y = rng.choice([-10000, -1000, -100, -17, -3, -2, -1, 0, 1, 2,
                        3, 17, 100, 1000, 10000, 2**53-1, 2**53,
                        2**53+2, 2**62, 2**63-1024])
        add("wide_integer", x, y)
    for _ in range(500):
        x = math.ldexp(1 + rng.random(), rng.randrange(-1074, 1024))
        y = rng.uniform(-30, 30)
        add("noninteger", x, y)
    for x in (0.71, 1.29):
        for neighbor in (math.nextafter(x, 0.), x, math.nextafter(x, math.inf)):
            for y in (-2**31-1, 2**31, -0.3, 0.3, 0.5):
                add("log_branch_threshold", neighbor, y)
    for x in (2., 0.5, 1e300, 1e-300, math.nextafter(1.,2.),
              math.nextafter(1.,0.)):
        for y in (-1e308, 1e308):
            add("huge_exponent_chain", x, y)
    for k in (510, 511, 512, 513, 537, 538):
        center=math.ldexp(1.0,k)
        for x in (math.nextafter(center,0.), center, math.nextafter(center,math.inf),
                  math.ldexp(1.2553091649534536,k)):
            for sign in (-1,1):
                for y in (-1, -2, -3):
                    add("reciprocal_restart_boundary", sign*x, y)
    for x in (math.nan, -math.nan, 0., -0., math.inf, -math.inf, 1., -1.):
        for y in (math.nan, -math.nan, math.inf, -math.inf, -0., 0., -0.5,
                  0.5, -1., 1., -3., 3., -2**31-1, -2**31, 2**31,
                  2**31+1, 2**53-1, 2**53):
            add("nonfinite_or_signed_zero", x, y)
    return cases

class ControlWord:
    """Tiny test-only Windows x86_64 reader/writer with explicit memory lifetime."""
    def __init__(self):
        self.kernel=ctypes.WinDLL("kernel32",use_last_error=True)
        self.kernel.VirtualAlloc.argtypes=[ctypes.c_void_p,ctypes.c_size_t,ctypes.c_ulong,ctypes.c_ulong]
        self.kernel.VirtualAlloc.restype=ctypes.c_void_p
        self.kernel.VirtualProtect.argtypes=[ctypes.c_void_p,ctypes.c_size_t,ctypes.c_ulong,ctypes.POINTER(ctypes.c_ulong)]
        self.kernel.VirtualFree.argtypes=[ctypes.c_void_p,ctypes.c_size_t,ctypes.c_ulong]
        self.kernel.FlushInstructionCache.argtypes=[ctypes.c_void_p,ctypes.c_void_p,ctypes.c_size_t]
        self.kernel.GetCurrentProcess.restype=ctypes.c_void_p
        self.address=self.kernel.VirtualAlloc(None,4096,0x3000,0x04)
        if not self.address:
            raise ctypes.WinError(ctypes.get_last_error())
        # fnstcw [rcx];ret;fldcw [rcx];ret. Read/write-only helpers, no math.
        code=b"\xd9\x39\xc3\xd9\x29\xc3"
        ctypes.memmove(self.address,code,len(code))
        previous=ctypes.c_ulong()
        if not self.kernel.VirtualProtect(self.address,4096,0x20,ctypes.byref(previous)):
            error=ctypes.get_last_error()
            self.close()
            raise ctypes.WinError(error)
        self.kernel.FlushInstructionCache(self.kernel.GetCurrentProcess(),self.address,len(code))
        signature=ctypes.CFUNCTYPE(None,ctypes.POINTER(ctypes.c_uint16))
        self.reader=signature(self.address)
        self.writer=signature(self.address+3)
    def get(self):
        result=ctypes.c_uint16()
        self.reader(ctypes.byref(result))
        return result.value
    def set(self,value):
        result=ctypes.c_uint16(value)
        self.writer(ctypes.byref(result))
    def close(self):
        if self.address:
            if not self.kernel.VirtualFree(self.address,0,0x8000):
                raise ctypes.WinError(ctypes.get_last_error())
            self.address=None

def primary(cases,control,mode):
    gmsh.initialize([],readConfigFiles=False)
    before=control.get()
    caller=(before&0xf0ff)|mode
    try:
        control.set(caller)
        for i,(_,x,y) in enumerate(cases):
            gmsh.parser.setNumber("x%d"%i,[x])
            gmsh.parser.setNumber("y%d"%i,[y])
        with tempfile.TemporaryDirectory(prefix="tessella_power_oracle_") as folder:
            script=Path(folder)/"inputs.geo"
            script.write_text("\n".join("v%d=x%d^y%d;"%(i,i,i)
                                         for i in range(len(cases))),encoding="utf8")
            gmsh.parser.parse(str(script))
        result=[bits(float(gmsh.parser.getNumber("v%d"%i)[0]))
                for i in range(len(cases))]
        assert control.get()==caller,"Gmsh changed caller control word"
        return result
    finally:
        control.set(before)
        gmsh.finalize()


def literal(x):
    if math.isnan(x):return "Sqrt(-1)" if math.copysign(1.,x)<0 else "-Sqrt(-1)"
    if math.isinf(x):return "Exp(10000)" if x>0 else "-Exp(10000)"
    return "%.17g"%x


def check_cli(cases,control,executable):
    groups=collections.defaultdict(list)
    for row in cases:groups[row[0]].append(row)
    selected=groups["signed_int32_cutoff"]+groups["noninteger"][:16]
    selected+=groups["log_branch_threshold"]+groups["huge_exponent_chain"]
    selected+=groups["reciprocal_restart_boundary"][:16]
    selected+=[row for row in groups["nonfinite_or_signed_zero"]
               if (row[1]==1 and math.isnan(row[2])) or
               (row[1] in (0.,-math.inf) and row[2] in
                (0.5,2147483649.,9007199254740991.,-2147483649.))]
    with tempfile.TemporaryDirectory(prefix="tessella_power_cli_") as folder:
        output=Path(folder)/"values.txt"
        source=['Printf(StrCat("GMSH_VERSION=",General.Version)) > "%s";'
                %output.as_posix()]
        source+=['Printf("%%.17g",(%s)^(%s)) >> "%s";'
                 %(literal(x),literal(y),output.as_posix()) for _,x,y in selected]
        script=Path(folder)/"inputs.geo"
        script.write_text("\n".join(source),encoding="utf8")
        subprocess.run([str(executable),str(script),"-parse_and_exit","-nopopup",
                        "-noenv"],capture_output=True,text=True,check=True,timeout=60)
        lines=output.read_text().splitlines()
    assert lines[0]=="GMSH_VERSION=4.15.2",lines[:3]
    actual=[bits(float(value)) for value in lines[1:]]
    expected=primary(selected,control,0x0300)
    assert len(actual)==len(selected)
    nan=lambda n:n&0x7ff0000000000000==0x7ff0000000000000 and n&0xfffffffffffff
    bad=[i for i,(a,b) in enumerate(zip(actual,expected))
         if a!=b and not(nan(a) and nan(b))]
    assert not bad,("CLI/PC64 DLL mismatches",bad[:10])
    print("standalone_cli_cases=%d mismatches=0 executable_sha256=%s"
          %(len(selected),sha(executable)))


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate",action="store_true")
    parser.add_argument("--gmsh-executable",type=Path)
    args=parser.parse_args()
    assert os.name=="nt" and platform.machine().lower() in ("amd64","x86_64")
    assert gmsh.__version__=="4.15.2",gmsh.__version__
    if args.gmsh_executable:
        assert args.gmsh_executable.is_file(),"native Gmsh executable does not exist"
        assert args.gmsh_executable.suffix.lower()==".exe",(
            "CLI oracle requires native gmsh.exe; Python gmsh.bat inherits hosted DLL precision")
    if not args.regenerate:
        metadata=json.loads((HERE/"corpus_metadata.json").read_text(encoding="utf8"))
        assert sha(gmsh.lib._name)==metadata["gmsh_library_sha256"],"primary DLL hash differs"
        assert sha(HERE/"corpus.tsv")==metadata["corpus_sha256"],"corpus hash differs"
        if args.gmsh_executable and metadata.get("gmsh_executable_sha256"):
            assert sha(args.gmsh_executable)==metadata["gmsh_executable_sha256"],(
                "primary executable hash differs")
    cases=case_catalog()
    control=ControlWord()
    try:
        original=control.get()
        columns=[primary(cases,control,mode) for mode in MODES]
        assert control.get()==original
        if args.gmsh_executable:check_cli(cases,control,args.gmsh_executable)
    finally:
        control.close()
    header="# group\tx_bits\ty_bits\t"+"\t".join("cw_%04x"%mode for mode in MODES)
    rows=[header]
    for i,(group,x,y) in enumerate(cases):
        rows.append("\t".join([group,"%016x"%bits(x),"%016x"%bits(y)]+
                              ["%016x"%column[i] for column in columns]))
    content=("\n".join(rows)+"\n").encode("ascii")
    corpus=HERE/"corpus.tsv"
    if args.regenerate:
        corpus.write_bytes(content)
        metadata={"gmsh_version":gmsh.__version__,"gmsh_library_sha256":sha(gmsh.lib._name),
                  "gmsh_executable_sha256":sha(args.gmsh_executable) if args.gmsh_executable else None,
                  "case_count":len(cases),"precision_bits":[53,64,24],
                  "control_word_mode_bits":MODES,"random_seed":221097,
                  "corpus_sha256":sha(corpus),"primary_script_sha256":sha(__file__),
                  "nan_comparison":"bitwise for DLL corpus; category after CLI text output",
                  "reference":"Gmsh parser x^y; Python assigns exact Float64 input bits",
                  "source_urls":[
                      "https://raw.githubusercontent.com/mingw-w64/mingw-w64/master/mingw-w64-crt/math/x86/pow.def.h",
                      "https://raw.githubusercontent.com/mingw-w64/mingw-w64/master/mingw-w64-crt/math/powi.def.h",
                      "https://raw.githubusercontent.com/mingw-w64/mingw-w64/master/mingw-w64-crt/math/x86/log2l.S",
                      "https://raw.githubusercontent.com/mingw-w64/mingw-w64/master/mingw-w64-crt/math/x86/exp2l.S"]}
        (HERE/"corpus_metadata.json").write_text(json.dumps(metadata,indent=2)+"\n",encoding="utf8")
    else:
        expected=corpus.read_bytes()
        assert content==expected,"pinned primary differs from committed corpus"
    print("primary_cases=%d modes=%d mismatches=0 corpus_sha256=%s"
          %(len(cases),len(MODES),sha(corpus)))

if __name__=="__main__":main()
