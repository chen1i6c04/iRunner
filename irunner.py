#!/usr/bin/env python3
import os
import sys
import shutil
import logging
import argparse
import xml.etree.cElementTree as ET
from utils import ApplicationError, run_cmd


def trimming(reads_1, reads_2, trim_1, trim_2, threads, length=36, min_qual=3):
    cmd = f"fastp -i {reads_1} -I {reads_2} -o {trim_1} -O {trim_2} --length_required {length} --cut_front {min_qual} "\
          f"--cut_tail {min_qual} --thread {threads} --detect_adapter_for_pe -j /dev/null -h /dev/null -t 1"
    run_cmd(cmd)


def dump_fastq(sra_file, outdir):
    run_cmd(f"fastq-dump --split-e --clip --outdir  {outdir} {sra_file}")


class SequenceReadArchive:
    def __init__(self, sra_run):
        self._check_input(sra_run)
        self._get_statistics()

    def _check_input(self, infile):
        with open(infile, 'rb') as handle:
            head = next(handle)[:8].decode()
        if head != 'NCBI.sra':
            raise TypeError(f"{infile} is not sra file.")
        else:
            self.run = infile

    def _get_statistics(self):
        stdout, stderr = run_cmd(f'sra-stat -xse 2 {self.run}')
        self._stat_tree = ET.fromstring(stdout)

    def count_bases(self, min_score=0):
        c = 0
        for quality in self._stat_tree.findall('*/Quality'):
            if int(quality.attrib['value']) >= min_score:
                c += int(quality.attrib['count'])
        return c

    @property
    def total_bases(self):
        return self.count_bases()

    @property
    def layout(self):
        return self._stat_tree.find('Statistics').attrib['nreads']

    @property
    def length(self):
        return tuple(map(lambda x: int(x.attrib['average']), self._stat_tree.findall('*//Read')))


def main():
    parser = argparse.ArgumentParser("NCBI SRA assembly pipeline.")
    parser.add_argument("sra", help="path of NCBI sra file")
    parser.add_argument("--outdir", required=True, help="Output folder")
    parser.add_argument("--tmpdir", default="/tmp", help="Directory of temp folder default: '/tmp'")
    parser.add_argument("--gsize", default='',
                        help="Estimated genome size(MB) eg. 3.2M, If blank will AUTODETECT. default: ''")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    parser.add_argument("--check", action='store_true', help="Check format of sra.")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    reads = os.path.join(args.outdir, 'READS.sra')
    os.symlink(args.sra, reads)
    logfile = os.path.join(args.outdir, 'runner.log')
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        format='[%(levelname)s] %(asctime)s %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    if args.check:
        try:
            sra = SequenceReadArchive(reads)
        except ApplicationError as exception:
            logging.debug(exception)
            sys.exit()
        if sra.layout != '2':
            logging.debug("Layout is not pair-end")
            sys.exit("Aborted")
        if sra.count_bases(30)/sra.count_bases()*100 < 80:
            logging.debug("q30 bases < 80%%")
            sys.exit("Aborted")
    fastq_dirname = os.path.join(args.outdir, 'fastq')
    logging.info("Dump fastq from SRA")
    try:
        dump_fastq(reads, fastq_dirname)
    except ApplicationError as exception:
        logging.debug(exception)
        sys.exit("Aborted!")
    raw_1 = os.path.join(fastq_dirname, 'READS_1.fastq')
    raw_2 = os.path.join(fastq_dirname, 'READS_2.fastq')
    paired_1 = os.path.join(fastq_dirname, 'R1.fastq.gz')
    paired_2 = os.path.join(fastq_dirname, 'R2.fastq.gz')
    logging.info('reads trimming')
    try:
        trimming(raw_1, raw_2, paired_1, paired_2, args.threads)
    except ApplicationError as exception:
        logging.debug(exception)
        sys.exit("Aborted!")
    logging.info("assembly with shovill")
    cmd = f"shovill --R1 {paired_1} --R2 {paired_2} --outdir {args.outdir} --depth 80 --tmpdir {args.tmpdir} " \
          f"--cpus {args.threads} --force --noreadcorr --nostitch"
    if args.gsize:
        cmd += f" --gsize {args.gsize}"
    try:
        run_cmd(cmd)
    except ApplicationError as exception:
        logging.debug(exception)
        sys.exit("Aborted!")
    shutil.rmtree(fastq_dirname)
    logging.info("Done")


if __name__ == '__main__':
    main()
