# parse_fastq.py

"""
Name:       parse_fastq.py
Author:     CAG
Version:    2.0
Date:       2026/03/24
"""

# %% Imports

import gzip
import regex
from tqdm import tqdm
from typing import Generator, Tuple

# %% Streaming functions

def stream_fastq(file_path: str) -> Generator[Tuple[str, str, str, str], None, None]:
    """
    Iteratively stream fq or fq.gz entries into memory as 4-line string.
    Shows an indeterminate tqdm progress bar (no pre-counting).
    """
    gz = file_path.endswith(".gz")

    if gz:
        # Buffer the underlying fileobj; gzip itself then reads bigger chunks
        raw = open(file_path, "rb", buffering=1024 * 1024)
        fh = gzip.open(raw, "rt", newline="")
    else:
        fh = open(file_path, "rt", buffering=1024 * 1024, newline="")

    with fh as fastq_file, tqdm(desc="Counting index/barcode pairs", unit="rec", ascii="-=") as pbar:
        it = iter(fastq_file)
        update_every = 1024
        n = 0

        while True:
            h = next(it, None)
            if h is None:
                break
            s = next(it)
            p = next(it)
            q = next(it)

            # cheaper than .strip(); we only want to drop newline
            yield h.rstrip("\n"), s.rstrip("\n"), p.rstrip("\n"), q.rstrip("\n")

            n += 1
            if n % update_every == 0:
                pbar.update(update_every)

        # flush remainder
        rem = n % update_every
        if rem:
            pbar.update(rem)


def process_fastq_stream(settings: dict):
    """
    Use stream_fastq() and parse_fastq_rec class to collect {idx:{bc:count}}
    """
    res_pfq = {}

    # Initialize search patterns
    parse_fastq_rec.init_patterns(settings)

    for _, seq, _, qual in stream_fastq(settings["sample_path"]):
        res = parse_fastq_rec(
            seq=seq,
            qual=qual,
            settings=settings,
        )

        # Ignore read if barcode is missing or low mean q
        if not hasattr(res, "bc_seq") or res.bc_qual < settings["mean_qual"]:
            continue

        # Process p5 - skip if missing or low mean q
        if not (hasattr(res, "p5_seq") and hasattr(res, "p5_qual")):
            continue

        p5_seq = res.p5_seq if res.p5_qual >= settings["mean_qual"] else None
        if p5_seq is None:
            continue

        # Optionally process p7, fill instead of skip:
        # Ns where q<mask_qual, or all Ns if mean q<min_qual, X if missing
        p7_seq = ""

        if settings["dualindex"]:
            if hasattr(res, "p7_seq") and hasattr(res, "p7_qual"):
                p7_seq = res.p7_seq if res.p7_qual >= settings["mean_qual"] else "N" * settings["tlen_p7"]
            else:
                p7_seq = "X" * settings["tlen_p7"]

        idx_seq = p5_seq + p7_seq

        if idx_seq not in res_pfq:
            res_pfq[idx_seq] = {}

        if res.bc_seq not in res_pfq[idx_seq]:
            res_pfq[idx_seq][res.bc_seq] = 0

        res_pfq[idx_seq][res.bc_seq] += 1

    return res_pfq


# %% Helpers

def mean_q(quals: list[int]) -> float:
    """
    Tiny hot-path helper: avoid NumPy overhead on short per-read quality slices.
    """
    return 0.0 if not quals else sum(quals) / len(quals)


# %% Class definition

class parse_fastq_rec(object):
    """
    Class to handle each fastq record as an instance while it streams.
    Return per-read idx (named if expected), bc, and assoc mean qualities.
    """

    __slots__ = (
        "seq",
        "qual",
        "p5_seq",
        "p5_qual",
        "bc_seq",
        "bc_qual",
        "p7_seq",
        "p7_qual",
    )

    # Class-level constants — built once, shared globally
    _RC_TABLE = None
    _PAT_P5 = None
    _PAT_P7 = None
    _PAT_BC = None

    @classmethod
    def init_patterns(cls, settings: dict):
        """
        Call once before streaming to precompile fuzzy patterns.
        This is much faster.
        """
        
        # Get allowed mismatches
        mm = settings["mismatches"]

        # compile regex patterns once; BESTMATCH yields the minimal-edit match
        cls._PAT_P5 = regex.compile(
            rf"(?:{regex.escape(settings['ref_p5'])}){{s<={mm}}}"   # FYI, e<={mm} & flags=regex.BESTMATCH takes somewhat longer but may be cleaner...
            )
        
        cls._PAT_BC = regex.compile(
            rf"(?:{regex.escape(settings['ref_bc'])}){{s<={mm}}}"
            )
        
        if settings.get("dualindex"):
            cls._PAT_P7 = regex.compile(
                rf"(?:{regex.escape(settings['ref_p7'])}){{s<={mm}}}"
                )

    def __init__(
        self, 
        seq: str, 
        qual: str, 
        settings: dict, 
        ):

        # Initialize sequence and quality
        self.seq = seq
        self.qual = self.Q33conv(qual)

        # (Opt) mask lowQ bases
        if settings["mask"]:
            self.seq = "".join(
                "N" if q < settings["mask_qual"] else b for b, q in zip(self.seq, self.qual)
                )

        # Process the fastq record according to settings
        self.process_sequence(settings)

    def process_sequence(self, settings):
        """
        Extract regions of interest and assoc. quality strings and update class attrs.
        Final attrs:
        - seq, qual
        - p5_seq, p5_qual
        - bc_seq, bc_qual
        - (Opt) p7_seq, p7_qual
        """
        # Extract p5 index
        p5 = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            pat=parse_fastq_rec._PAT_P5,
            targ_dir=settings["tdir_p5"],
            targ_len=settings["tlen_p5"],
            )

        # If p5 found, set attrs and continue
        if p5[0]:
            self.p5_seq, self.p5_qual = p5[0], mean_q(p5[1])

            # Reverse complement if called
            if settings["rdir"] == "rev":
                self.seq = self.revcomp(self.seq)
                self.qual = self.qual[::-1]

            # Process as single or dual index
            if not settings["dualindex"]:
                self.process_single_index(settings)
            else:
                self.process_dual_index(settings)

    def process_single_index(self, settings):
        """
        Obtain barcode based on target and length
        """
        # Extract barcode for single index
        bc = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            pat=parse_fastq_rec._PAT_BC,
            targ_dir=settings["tdir_bc"],
            targ_len=settings["tlen_bc"],
            )

        if bc[0]:
            self.bc_seq, self.bc_qual = bc[0], mean_q(bc[1])

    def process_dual_index(self, settings):
        """
        Obtain p7, and if present obtain barcode based on two targets to 
        facilitate optional repair of short barcode, given tlen_bc and a fill sequence.
        """
        # Get di_profile
        base_profile = settings.get("base_profile")

        # Extract p7 index
        p7 = self.match_extract_single_targ(
            seq=self.seq,
            qual=self.qual,
            pat=parse_fastq_rec._PAT_P7,
            targ_dir=settings["tdir_p7"],
            targ_len=settings["tlen_p7"],
            )

        # Exit loop if no p7
        if not p7[0]:
            return
        
        # If found, set attrs and continue
        self.p7_seq, self.p7_qual = p7[0], mean_q(p7[1])

        # Finish based on di_profile
        # Extract barcode for dual index:
        
        if base_profile == "X": 
            # `X` has a long amplicon, just extract downstream from (pre-compiled) ref_bc. Also...
            bc = self.match_extract_single_targ( 
                seq=self.seq,
                qual=self.qual,
                pat=parse_fastq_rec._PAT_BC,
                targ_dir=settings["tdir_bc"],
                targ_len=settings["tlen_bc"],
                )
            # The X/INT p7 seq is revcomp, but bc is forward. I don't like this fix, but it works.
            self.p7_seq = self.revcomp(self.p7_seq)
            
        else:
            bc = self.match_extract_dual_targ( #`M`s are flanked by PAT_P7 and PAT_BC, cut between
                seq=self.seq,
                qual=self.qual,
                pat1=parse_fastq_rec._PAT_P7,
                pat2=parse_fastq_rec._PAT_BC,
                targ_len=settings["tlen_bc"],
                fill_seq=settings["fill_seq"],
                )

        if bc[0]:
            self.bc_seq, self.bc_qual = bc[0], mean_q(bc[1])

    @staticmethod
    def Q33conv(qual: str):
        """
        Convert a string of Q33 ascii characters into a list of PHRED scores.
        """
        return [ord(q) - 33 for q in qual]

    @staticmethod
    def revcomp(seq: str):
        """
        Returns the reverse complement of a DNA sequence.
        Lazily initialize the RC table if needed.
        """
        if parse_fastq_rec._RC_TABLE is None:
            parse_fastq_rec._RC_TABLE = str.maketrans(
                {"A":"T","C":"G","G":"C","T":"A","N":"N","X":"X"}
            )
        return seq.translate(parse_fastq_rec._RC_TABLE)[::-1]

    @staticmethod
    def match_extract_single_targ(
        seq: str, 
        qual: list[int], 
        pat: "regex.Pattern", 
        targ_dir: str, 
        targ_len: int,
        ):

        """
        Extract sequence of given length and direction by single target.
        """
        # locate match
        m = pat.search(seq)

        if m:
            # get positions
            start, end = m.start(), m.end()

            # return appropriate slice
            if targ_dir == "upstream":
                return [seq[start - targ_len : start], qual[start - targ_len : start]]
            elif targ_dir == "downstream":
                return [seq[end : end + targ_len], qual[end : end + targ_len]]

            # or tell user we need a direction
            else:
                raise ValueError("targ_dir param requires 'upstream' or 'downstream'")

        # if no match, return empty
        return ["", []]

    @staticmethod
    def match_extract_dual_targ(
        seq: str,
        qual: list[int],
        pat1: "regex.Pattern",
        pat2: "regex.Pattern",
        targ_len: int,
        fill_seq: str | None,
        ):
        
        """
        Extract sequence between ref1 and ref2.
        Pad front of xseq to targ_len with 3' of fill_seq, if fill_seq provided
        """
        # locate matches
        m1 = pat1.search(seq)
        m2 = pat2.search(seq)

        if m1 and m2:
            start, end = m1.end(), m2.start()

            # if start is before or at end (this catches empty short bcs!)
            if start <= end:
                # get slices
                xseq, xqual = seq[start:end], qual[start:end]

                # pad 5' bc_seq to desired length with 3' of fill_seq
                # (if fill_seq is provided in stock settings, and bc is short)
                if fill_seq and len(xseq) < targ_len:
                    diff = targ_len - len(xseq)
                    xseq = fill_seq[-diff:] + xseq
                    xqual = [30] * diff + xqual

                return [xseq, xqual]

        return ["", []]


#%% Versions
"""
2.0
- pipeline v4 overhaul, all fastq processing moved here 
1.5 
- added logic for dual-index X/INT
1.4 
- pre-compiled regex patterns, new parameter passing method via settings
"""
