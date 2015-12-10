from __future__ import print_function
import collections
import os.path
import re
import shutil

from bs4 import BeautifulSoup
import requests

LOCI = collections.OrderedDict()
LOCI['A'] = {'Reference': '01:01:01:01', 'Type': 'Genomic'}
LOCI['B'] = {'Reference': '07:02:01', 'Type': 'Genomic'}
LOCI['C'] = {'Reference': '01:02:01', 'Type': 'Genomic'}
LOCI['DPA1'] = {'Reference': '01:03:01:01', 'Type': 'Genomic'}
LOCI['DPB1'] = {'Reference': '01:01:01', 'Type': 'Genomic'}
LOCI['DQA1'] = {'Reference': '01:01:01', 'Type': 'Genomic'}
LOCI['DQB1'] = {'Reference': '05:01:01:01', 'Type': 'Genomic'}
LOCI['DRB1'] = {'Reference': '01:01:01', 'Type': 'Genomic'}
LOCI['DRB3'] = {'Reference': '01:01:01', 'Type': 'Genomic'}
LOCI['DRB4'] = {'Reference': '01:01:01', 'Type': 'Genomic'}
LOCI['DRB5'] = {'Reference': '01:01:01', 'Type': 'CDS'}


def locus_path(locus):
    return 'downloaded/{}.html'.format(locus)


NOT_WHITESPACE = re.compile(r'\S', re.UNICODE)


class OrderedDictOfLists(collections.OrderedDict):
    def extend(self, key, elements):
        if key not in self:
            self[key] = []
        self[key].extend(elements)


for locus in LOCI:
    if os.path.exists(locus_path(locus)):
        print('Already downloaded', locus)
        continue
    print('Downloading', locus)
    # Somehow, this webserver checks the order the post data is in, so use an
    # OrderedDict to force iteration to occur in the same order as insertion
    post_data = collections.OrderedDict()
    post_data['gene'] = locus
    post_data['Type'] = LOCI[locus]['Type']
    post_data['Reference'] = LOCI[locus]['Reference']
    post_data['Sequences'] = ''
    post_data['Display'] = 'Show All Bases'
    post_data['Formatting'] = 10
    post_data['Omit'] = 'N'
    post_data['Printing'] = 'P'
    post_data['submit'] = 'Align Sequences Now'
    locus_page = requests.post(
        'http://www.ebi.ac.uk/cgi-bin/ipd/imgt/hla/align.cgi', data=post_data,
        stream=True)
    with open('tmp', 'wb') as f:
        for chunk in locus_page.iter_content(chunk_size=1024*1024):
            f.write(chunk)
            print('.', end='')
    print()
    shutil.move('tmp', locus_path(locus))
    print('Downloaded', locus)

for locus in LOCI:
    print('Processing', locus)
    tree = BeautifulSoup(open(locus_path(locus), 'rb'), 'html.parser')
    header_line = None
    next_line_is_header = False
    # Store the sequences in the same order they appear in the document
    rows = OrderedDictOfLists()
    for element in tree.find('pre').find_all(string=True):
        line = unicode(element.string)
        if len(line) == 1:
            # Every block of lines starts with a single-space line, so the next
            # line is the header
            next_line_is_header = True
            continue
        if next_line_is_header:
            # The header line has to be processed differently depending on how
            # lines in the body are processed, so don't split it right now
            header_line = line
            next_line_is_header = False
            continue
        columns = line.split()
        # For the first body line after the header line, use the position of the
        # splits in that body line to figure out where the columns in the header
        # line are
        if header_line is not None:
            column_end = 0
            header = []
            header_push = 0
            for column in columns:
                # Start looking for the next column at the end of the previous
                # column
                column_start = line.index(column, column_end)
                column_end = column_start + len(column) - 1
                header_start = header_push + column_start
                header_end = header_push + column_end
                # Normally, the column header ends at "column_end", but there
                # are columns that aren't wide enough for the column header, so
                # the rest of the column headers get pushed forwards.
                while (header_end + 1 < len(header_line) - 1) and (
                        NOT_WHITESPACE.match(header_line[header_end + 1])):
                    header_end += 1
                header_push += (header_end - column_end)
                header.append(
                    header_line[header_start:header_end + 1].strip())
            rows.extend(header[0], header[1:])
            header_line = None
        rows.extend(columns[0], columns[1:])
    with open('{}-split.csv'.format(locus), 'w') as f:
        for row in rows:
            f.write('{},{}\n'.format(row, ','.join(rows[row])))
    with open('{}-combined.csv'.format(locus), 'w') as f:
        first_row = True
        for row in rows:
            # Skip the first row since column headers don't make sense without
            # splitting data into columns
            if first_row:
                first_row = False
                continue
            f.write('{},{}\n'.format(row, ''.join(rows[row])))
