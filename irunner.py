#!/usr/bin/env python3
import os
import re
import sys
import shutil
import logging
import argparse
import xml.etree.cElementTree as ET
from subprocess import run, PIPE, CalledProcessError


def run_cmd(cmd):
    p = run(cmd, stdout=PIPE, stderr=PIPE, shell=True, check=True)
    return p


def trimming(reads_1, reads_2, trim_1, trim_2, threads):
    cmd = f"fastp -i {reads_1} -I {reads_2} -o {trim_1} -O {trim_2} --length_required 36 --cut_front 3 --cut_tail 3 " \
          f"--thread {threads} --detect_adapter_for_pe -j /dev/null -h /dev/null"
    run_cmd(cmd)


def download_and_dump_fastq(accession, outdir):
    run_cmd(f"fastq-dump --split-e --outdir {outdir} {accession}")


class SequenceReadArchive:
    def __init__(self, accession):
        self._parse_accession(accession)
        self._get_statistics()

    def _parse_accession(self, accession):
        if re.fullmatch('^[DES]RR[0-9]+$', accession):
            self._run_accession = accession
        else:
            raise ValueError(f"{accession} is not run accession")

    def _get_statistics(self):
        process = run_cmd(f'sra-stat -xse 2 {self._run_accession}')
        self._stat_tree = ET.fromstring(process.stdout.decode())

    def count(self, min_score=0):
        c = 0
        for quality in self._stat_tree.findall('*/Quality'):
            if int(quality.attrib['value']) >= min_score:
                c += int(quality.attrib['count'])
        return c

    @property
    def total_bases(self):
        return self.count()

    @property
    def layout(self):
        return self._stat_tree.find('Statistics').attrib['nreads']

    def high_qulity_bases_percent(self, qscore=30):
        return self.count(qscore) / self.count() * 100

    @property
    def length(self):
        return tuple(map(lambda x: int(x.attrib['average']), self._stat_tree.findall('*//Read')))


def main():
    parser = argparse.ArgumentParser("NCBI SRA assembly pipeline.")
    parser.add_argument("accession", help="NCBI Run accession")
    parser.add_argument("--outdir", required=True, help="Output folder")
    parser.add_argument("--tmpdir", default="/tmp", help="Directory of temp folder default: '/tmp'")
    parser.add_argument("--gsize", default='',
                        help="Estimated genome size(MB) eg. 3.2M, If blank will AUTODETECT. default: ''")
    parser.add_argument("--threads", default=8, type=int, help="Number of threads to use. default: 8")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    logfile = os.path.join(args.outdir, 'runner.log')
    logging.basicConfig(
        filename=logfile,
        level=logging.DEBUG,
        format='[%(levelname)s] %(asctime)s %(message)s',
        datefmt="%Y-%m-%d %H:%M:%S"
    )

    try:
        sra = SequenceReadArchive(args.accession)
    except Exception as e:
        logging.debug(e)
        sys.exit()
    if sra.layout != '2':
        logging.debug("Layout is not pair-end")
        sys.exit("Aborted")
    if sra.high_qulity_bases_percent() < 80:
        logging.debug("q30 bases < 80%%")
        sys.exit("Aborted")
    fastq_dirname = os.path.join(args.outdir, 'fastq')
    logging.info("Dump fastq from SRA")
    try:
        download_and_dump_fastq(args.accession, fastq_dirname)
    except CalledProcessError:
        logging.debug("Dump fastq fail")
        sys.exit("Aborted!")
    reads_1 = os.path.join(fastq_dirname, f'{args.accession}_1.fastq')
    reads_2 = os.path.join(fastq_dirname, f'{args.accession}_2.fastq')
    trim_reads_1 = os.path.join(fastq_dirname, 'R1.fastq')
    trim_reads_2 = os.path.join(fastq_dirname, 'R2.fastq')
    logging.info('reads trimming')
    try:
        trimming(reads_1, reads_2, trim_reads_1, trim_reads_2, args.threads)
    except CalledProcessError:
        logging.debug("reads trimming fail")
        sys.exit("Aborted!")
    logging.info("assembly with shovill")
    cmd = f"shovill --R1 {trim_reads_1} --R2 {trim_reads_2} --outdir {args.outdir} --depth 80 --tmpdir {args.tmpdir} " \
          f"--cpus {args.threads} --force --noreadcorr --nostitch"
    if args.gsize:
        cmd += f" --gsize {args.gsize}"
    try:
        run_cmd(cmd)
    except CalledProcessError:
        logging.debug("assembly fail")
        sys.exit("Aborted!")
    shutil.rmtree(fastq_dirname)
    logging.info("Done")


if __name__ == '__main__':
    main()
