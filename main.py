from argparse import ArgumentParser
from itertools import product
from collections import namedtuple, defaultdict
import random
import logging

# Helpers
_COMMENT = ';'
_RC_PREF = '_'
_NUCS = 'a c g t'.upper().split()
_RCs =  't g c a'.upper().split()
_RCs =  dict(zip(_NUCS, _RCs))

def rand_kmer(k):
    return ''.join(random.choice(_NUCS) for i in range(k))

def rc(seq):
    return ''.join(_RCs[c] for c in seq)[::-1]

def to_cannonical(km):
    kw = rc(km)
    return km if km < kw else kw

def exit_failed():
    logging.error('Failed to find compatible contig sequences.')
    logging.error('Try again, or try with larger k')
    exit(1)

###

Contig = namedtuple('Contig', 'id o')

class RefSeqFromTilings:
    def __init__(self, tilings, k):
        '''
        A tiling to reference sequence generator.
        Tilings represented as list of (Contig ID: str, Orientation: bool) pairs
        '''
        self.k = k
        self.tilings = tilings

        self.contigs = {}
        self.ref_seqs = []

        self.kmers = set()

        self._heads = defaultdict(lambda: None) # prefixes
        self._tails = defaultdict(lambda: None) # suffixes

        for tiling in self.tilings:
            for c in tiling:
                self.contigs[c.id] = None

    # def rc_edge(c):
    #     return (Contig(c[1].id, not c[0].o), Contig(c[0].id, not c[1].o))
    @classmethod
    def from_file(cls, fp, k):
        '''
        Generate from tilings file.
        Each reference on a new line. Comment lines starting with ';'.
        Tilings represented space seperated sequence of <o><ctg_id> pairs
        o = '_' if rc and '' otherwise
        '''
        tilings = []
        with open(fp) as f:
            for line in f:
                if line.startswith(_COMMENT): continue
                ref_contigs = line.split()
                ref_contigs = [Contig(c[1:], False) if c.startswith(_RC_PREF)
                               else Contig(c, True) for c in ref_contigs]
                tilings.append(ref_contigs)

        obj = cls(tilings, k=k)
        return obj

    def generate(self):
        '''
        Generate contigs and reference sequences from tiling
        1. Iterate over tilings in sequence inserting random k-1 overlaps
        2. Iterate over contigs with (k-1) suffix and prefixes and insert
           nucleotides to guarantee that no contigs share kmers
        3. Read out contig tilings to generate tilings
        '''
        self.generate_contigs()
        for tiling in self.tilings:
            seqs = []
            for ctg in tiling:
                c_seq = self.contigs[ctg.id]
                if not ctg.o:
                    c_seq = rc(c_seq)
                c_seq = c_seq[self.k - 1:] if seqs else c_seq
                seqs.append(c_seq)
            self.ref_seqs.append(''.join(seqs))

    def generate_contig_overlaps(self):
        '''
        Generate contig (k-1) suffixes and prefixes.
        These must be consistent with contig-contig k-1 as defined by reference tilings.
        1. Iterate over adjacent pairs on tilings. To determine (k-1) overlap.
        2a. If neither contig in contig pair has an overlap defined by suff/prefix,
            then generate a random (k-1) mer for the overlap, and fill in
            'terminal' k-1 bases (accounting for orientation)
        2b. If one of the contig pairs has (k-1) terminal bases participating in then
            overlap defined, then propagate those bases to the other contig
        2c. If both contig pairs has (k-1) terminla bases participating in the overlap
            defined, then assert that those (k-1) bases are equal (up to orientation)
        '''

        for tiling in self.tilings:
            preds = tiling[:-1]
            succs = tiling[1:]
            for p, s in zip(preds, succs):
                p_overlap = self._tails[p.id] if p.o else rc(self._heads[p.id])
                s_overlap = self._heads[s.id] if s.o else rc(self._tails[s.id])
                if p_overlap is None:
                    if s_overlap is None:
                        overlap = rand_kmer(self.k - 1)
                    else:
                        overlap = s_overlap
                else:
                    if s_overlap is None:
                        overlap = p_overlap
                    else:
                        assert(p_overlap == s_overlap)
                        overlap = p_overlap

                if p.o:
                    self._tails[p.id] = overlap
                else:
                    self._heads[p.id] = rc(overlap)

                if s.o:
                    self._heads[s.id] = overlap
                else:
                    self._tails[s.id] = rc(overlap)

    def _is_valid_new_ctg(self, ctg):
        '''
        Determines if given contig can be a valid new contig given
        all other contigs.
        '''
        if len(ctg) < self.k: return False
        ctg_kmers = set()
        end_pos = len(ctg) - self.k + 1
        for pos in range(end_pos):
            km = ctg[pos: pos + self.k]
            km = to_cannonical(km)
            if (km in ctg_kmers) or (km in self.kmers):
                return False
            ctg_kmers.add(km)
        return True

    def generate_contigs(self):

        def insert_ctg(ctg_id, ctg_seq):
            self.contigs[c] = ctg
            insert_ctg_kmers(ctg)

        def insert_ctg_kmers(ctg):
            end_pos = len(ctg) - self.k + 1
            for pos in range(end_pos):
                km = ctg[pos: pos + self.k]
                self.kmers.add(to_cannonical(km))

        def gen_new_contig(head, tail):
            # Try (lexicographically) to insert bases
            # to find new valid contig
            for l in range(0, self.k):
                for ins in product(_NUCS, repeat=l):
                    ins = ''.join(ins)
                    ctg = head + ins + tail
                    if self._is_valid_new_ctg(ctg):
                        return ctg
            exit_failed()

        # Generate terminal k-1 bases for contigs with neighbors on dbG
        self.generate_contig_overlaps()
        # Generate contig sequences one at a time
        for c in self.contigs.keys():
            head = self._heads[c]
            tail= self._tails[c]
            head = head if head else ''
            tail = tail if tail else ''
            ctg = gen_new_contig(head, tail)
            insert_ctg(c, ctg)

    def debug(self):
        tot_len = sum(len(ref) for ref in self.ref_seqs)
        logging.debug('Total reference length: {}'.format(tot_len))
        logging.debug('Generated compatible contigs with {} {}-mers'
                      .format(len(self.kmers), self.k))
        for c in sorted(self.contigs.keys()):
            logging.debug('\t{}: {}'.format(c, self.contigs[c]))

    def fasta(self):
        lines = []
        for i, s in enumerate(self.ref_seqs):
            lines.append('>{}'.format(i))
            lines.append(s.upper())
        return '\n'.join(lines) + '\n'

    def contigs_yaml(self):
        lines = []
        for c, seq in self.contigs.items():
            lines.append('{}: {}'.format(c, self.contigs[c]))
        return '\n'.join(lines) + '\n'

    def print_fasta(self):
        print(self.fasta())

    def print_contigs(self):
        print(self.contigs_yaml())

    def to_fasta_file(self, fp):
        with open(fp, 'w') as f:
            f.write(self.fasta())

    def to_contigs_file(self, fp):
        with open(fp, 'w') as f:
            f.write(self.contigs_yaml())

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-i', '--tilings', required=True)
    parser.add_argument('-v', '--verbose', required=False, default=1, type=int)
    parser.add_argument('-k', '--k', required=False, default=5, type=int)
    parser.add_argument('-o', '--out_dir', required=False, type=str)
    parser.add_argument('-n', '--name', required=False, type=str)
    parser.add_argument('--seed', required=False, default=2022, type=int)
    parser.add_argument('--print', required=False, action='store_true')
    parser.add_argument('--print_contigs', required=False, action='store_true')
    return parser.parse_args()

def main():
    args = parse_args()
    random.seed(args.seed)
    logging.basicConfig(level=args.verbose)

    refseqs = RefSeqFromTilings.from_file(args.tilings, args.k)
    refseqs.generate()
    refseqs.debug()
    refseqs.print_fasta()

    if args.out_dir:
        import os
        if not os.path.isdir(args.out_dir):
            os.mkdir(args.out_dir)

        name = args.name if args.name else args.out_dir
        pref = os.path.join(args.out_dir, name)
        refseqs.to_fasta_file(pref + '.fasta')
        refseqs.to_contigs_file(pref + '.ctgs')

if __name__ == "__main__":
    main()
