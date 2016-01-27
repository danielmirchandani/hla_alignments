#!/usr/bin/python
from __future__ import print_function
import collections
import os.path
import re
import shutil
import sys

import bs4
import requests

_LOCI = {
    'A': {'Reference': '01:01:01:01', 'Type': 'Genomic'},
    'B': {'Reference': '07:02:01', 'Type': 'Genomic'},
    'C': {'Reference': '01:02:01', 'Type': 'Genomic'},
    'DPA1': {'Reference': '01:03:01:01', 'Type': 'Genomic'},
    'DPB1': {'Reference': '01:01:01', 'Type': 'Genomic'},
    'DQA1': {'Reference': '01:01:01', 'Type': 'Genomic'},
    'DQB1': {'Reference': '05:01:01:01', 'Type': 'Genomic'},
    'DRB1': {'Reference': '01:01:01', 'Type': 'Genomic'},
    'DRB3': {'Reference': '01:01:01', 'Type': 'Genomic'},
    'DRB4': {'Reference': '01:01:01', 'Type': 'Genomic'},
    'DRB5': {'Reference': '01:01:01', 'Type': 'CDS'},
}

_LOCI_URL = 'http://www.ebi.ac.uk/cgi-bin/ipd/imgt/hla/align.cgi'

_NOT_WHITESPACE = re.compile(r'\S', re.UNICODE)

_OUTPUT_DIRECTORY = 'hla-alignments'


def _download_locus(locus, output_path):
    if os.path.exists(output_path):
        print('Already downloaded', locus)
        return
    print('Downloading', locus)
    # This webserver uses POST data by index, not by name; use an OrderedDict
    # so POST data is sent in this specific order.
    post_data = collections.OrderedDict()
    post_data['gene'] = locus
    post_data['Type'] = _LOCI[locus]['Type']
    post_data['Reference'] = _LOCI[locus]['Reference']
    post_data['Sequences'] = ''
    post_data['Display'] = 'Show All Bases'
    post_data['Formatting'] = 10
    post_data['Omit'] = 'N'
    post_data['Printing'] = 'P'
    post_data['submit'] = 'Align Sequences Now'
    locus_page = requests.post(_LOCI_URL, data=post_data, stream=True)
    with open('tmp', 'wb') as temp_file:
        for chunk in locus_page.iter_content(chunk_size=1024*1024):
            temp_file.write(chunk)
            print('.', end='', file=sys.stdout)
            sys.stdout.flush()
    print()
    shutil.move('tmp', output_path)
    print('Downloaded', locus)


def _get_lines_from_download(input_path):
    soup = bs4.BeautifulSoup(open(input_path, 'rb'), 'html.parser')
    # The text starts after the first <br>, so ignore any text before it
    current_br = soup.find('pre').br
    line = ''
    # BeautifulSoup treats <br> as enclosing all further content until the
    # tag enclosing the <br> is closed. For example, for "<p>a<br>b<br>c</p>",
    # BeautifulSoup creates:
    # Tag('p', ['a', Tag('br', ['b', Tag('br', ['c'])])])
    while current_br:
        next_br = None
        for element in current_br.children:
            if isinstance(element, bs4.element.Tag):
                if element.name == 'br':
                    # A <br> element means the end of the current line
                    yield line
                    line = ''
                    next_br = element
                elif element.name == 'span':
                    # Treat <span> elements as continuing the line
                    line += unicode(element.string)
                else:
                    print('Found tag with unknown name', element.name,
                          file=sys.stderr)
            elif isinstance(element, bs4.element.NavigableString):
                line += unicode(element.string)
            else:
                print('Found unknown type', type(element), file=sys.stderr)
        current_br = next_br
    yield line


def _process_locus(locus, input_path, output_path):
    print('Processing', locus)
    next_line_is_header = False
    rows = collections.defaultdict(lambda: [])
    for line in _get_lines_from_download(input_path):
        if len(line) == 1:
            # Every block of lines starts with a single-space line, so the next
            # line is the header
            next_line_is_header = True
            continue
        if next_line_is_header:
            # The header line is useless for this processing, so skip it
            next_line_is_header = False
            continue
        columns = line.split()
        rows[columns[0]].extend(columns[1:])
    with open(output_path, 'w') as combined_file:
        for row in sorted(rows):
            combined_file.write('{},{}\n'.format(row, ''.join(rows[row])))


def main():
    if not os.path.exists(_OUTPUT_DIRECTORY):
        os.mkdir(_OUTPUT_DIRECTORY)
    for locus in sorted(_LOCI):
        download_path = '{}/{}.html'.format(_OUTPUT_DIRECTORY, locus)
        _download_locus(locus, download_path)
        _process_locus(locus, download_path,
                       '{}/{}.csv'.format(_OUTPUT_DIRECTORY, locus))

if __name__ == '__main__':
    main()
