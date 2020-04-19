#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by huangzhibo on 2020/4/19
import gzip
import os
import sys
from argparse import ArgumentParser
import pysam


class AllsiteVcf:
    def __init__(self, vcf):
        self.input_vcf = vcf
        self.vcf_reader = pysam.VariantFile(self.input_vcf)
        # if not os.path.exists(self.vcf_reader.index):
        #     raise FileExistsError("No vcf index: {}".format(self.input_vcf))

    def get_genotype(self, chromosome: str, position: int):
        gt = '--'
        rec_it = self.vcf_reader.fetch(chromosome, position - 1, position)
        for rec in rec_it:
            # print(rec)
            # pysam.VariantRecord(rec).pos == position
            if rec.pos == position:
                first, second = rec.samples.get(0)['GT']
                if first is None:   # 跳过nocall位点 ./.
                    continue
                if rec.filter.keys() and rec.filter.keys()[0] != '.' and rec.filter.keys()[0] != 'PASS':
                    continue
                if rec.samples.get(0).get('GQ') and rec.samples.get(0).get('GQ') < 20:
                    continue
                if rec.samples.get(0).get('DP') and rec.samples.get(0).get('DP') < 10:
                    continue
                a1 = rec.alleles[first]
                a2 = rec.alleles[second]
                if len(rec.ref) == 1:
                    if len(a1) == 1:
                        a1_type = a1
                    else:
                        a1_type = 'I'
                    if len(a2) == 1:
                        gt = a1_type+a2
                    else:
                        gt = a1_type+'I'
                else:
                    a1_type = self.get_indel_alelle_type(rec.ref, a1)
                    gt = a1_type + self.get_indel_alelle_type(rec.ref, a2)
        return gt

    @staticmethod
    def get_indel_alelle_type(ref, alt):
        """
            CTTT    C,TTTT
            TGCGC   T,TGTGCGCGCGC
        :param ref: ref
        :param alt: alt
        :return: genotype
        """
        if len(ref) - len(alt) > 0:
            alele_type = 'D'
        elif len(ref) - len(alt) < 0:
            alele_type = 'I'
        else:
            alele_type = alt[0]
        return alele_type


def main():
    parser = ArgumentParser(description='Convert allsite VCF to wegene/23andMe raw data format')
    parser.add_argument("-i", "--input", dest="input", required=True,
                        help='allsite.vcf.gz with index [requested]')
    parser.add_argument("-b", "--blank", dest="blank", required=True,
                        help='wegene/23andMe format file without GT [requested]')
    parser.add_argument("-o", "--output", dest="output", default='out.txt', help='output wegene/23andMe format file')

    if len(sys.argv) == 1:
        parser.print_help()
        exit(0)
    args = parser.parse_args()
    if not os.path.exists(args.input):
        raise Exception('-i: The path of the file for input does not exist.')
    if not os.path.exists(args.blank):
        raise Exception('-b: The path of the file for blank does not exist.')

    avcf = AllsiteVcf(args.input)

    if os.path.splitext(args.blank)[1] == '.gz':
        bl = gzip.open(args.blank,'rt')
    else:
        bl = open(args.blank,'r')

    with open(args.output, 'w') as output_file:
        for line in bl:
            line = line.rstrip()
            output_file.write(line)

            if line.startswith("#") or not line:
                output_file.write('\n')
                continue

            identifier, chromosome, position = line.split("\t")[:3]
            gt = avcf.get_genotype(chromosome, int(position))
            output_file.write('\t')
            output_file.write(gt)
            output_file.write('\n')

    bl.close()


if __name__ == "__main__":
    sys.exit(main())
