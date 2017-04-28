#!/usr/bin/env python
"""
#
#------------------------------------------------------------------------------
#                         University of Minnesota
#         Copyright 2013, Regents of the University of Minnesota
#------------------------------------------------------------------------------
# Author:
#
#  James E Johnson
#
#------------------------------------------------------------------------------
"""


"""
This tool takes a SnpEff VCF file and an Ensembl pep.all.fa file ( e.g. Homo_sapiens.GRCh37.73.pep.all.fa )
It outputs a peptide fasta file with the variant peptide sequence that result from NON_SYNONYMOUS_CODING effects 

"""

import sys,re,os.path
import tempfile
import optparse
from optparse import OptionParser
import logging

## dictionary for Amino Acid Abbreviations
aa_abbrev_dict = dict()
aa_abbrev_dict['Phe'] = 'F'
aa_abbrev_dict['Leu'] = 'L'
aa_abbrev_dict['Ser'] = 'S'
aa_abbrev_dict['Tyr'] = 'Y'
aa_abbrev_dict['Cys'] = 'C'
aa_abbrev_dict['Trp'] = 'W'
aa_abbrev_dict['Pro'] = 'P'
aa_abbrev_dict['His'] = 'H'
aa_abbrev_dict['Gln'] = 'Q'
aa_abbrev_dict['Arg'] = 'R'
aa_abbrev_dict['Ile'] = 'I'
aa_abbrev_dict['Met'] = 'M'
aa_abbrev_dict['Thr'] = 'T'
aa_abbrev_dict['Asn'] = 'N'
aa_abbrev_dict['Lys'] = 'K'
aa_abbrev_dict['Val'] = 'V'
aa_abbrev_dict['Ala'] = 'A'
aa_abbrev_dict['Asp'] = 'D'
aa_abbrev_dict['Glu'] = 'E'
aa_abbrev_dict['Gly'] = 'G'

##  Get the peptide ID and sequence a given ID 
def get_sequence(id,seq_file):
  fh = open(seq_file, 'r')
  try:
    for (ln,line) in enumerate(fh):
      if line.find(id) >= 0:
        fields = line.split('\t')
        return ( ' '.join(fields[0:-1]),fields[-1].rstrip() if fields and len(fields) > 0 else None )
  except Exception, e:
    print >> sys.stderr, "failed: %s" % e
  finally:
    fh.close()

def fasta_to_tabular(fasta_file,tabular_file):
  inFile = open(fasta_file,'r')
  outFile = open(tabular_file,'w') 
  for i, line in enumerate( inFile ):
    line = line.rstrip( '\r\n' )
    if not line or line.startswith( '#' ):
      continue
    if line.startswith( '>' ):
      #Don't want any existing tabs to trigger extra columns:
      line = line.replace('\t', ' ')
      if i > 0:
        outFile.write('\n')
      outFile.write(line[1:])
      outFile.write('\t')
    else:
      outFile.write(line)
  if i > 0:
    outFile.write('\n')
  if inFile:
    inFile.close()
  if outFile:
    outFile.close()

def __main__():
  #Parse Command Line
  parser = optparse.OptionParser()
  parser.add_option( '-i', '--input', dest='input', help='The input snpeff vcf file with HGVS annotations (else read from stdin)' )
  parser.add_option( '-o', '--output', dest='output', help='The output fasta (else write to stdout)' )
  parser.add_option( '-p', '--protein_fasta', dest='protein_fasta', default=None, help='The Esembl protein fasta in tabular format' )
  parser.add_option( '-l', '--leading_aa_num', dest='leading_aa_num', type='int', default=None, help='leading number of AAs to output' )
  parser.add_option( '-t', '--trailing_aa_num', dest='trailing_aa_num', type='int', default=None, help='trailing number of AAs to output' )
  parser.add_option( '-d', '--debug', dest='debug', action='store_true', default=False, help='Turn on wrapper debugging to stdout'  )
  (options, args) = parser.parse_args()

  # need protein_fasta file
  fastaFile = options.protein_fasta
  if options.protein_fasta == None:
    print >> sys.stderr, "Ensembl protein_fasta tabular file required"
    exit(4)
  else:
    # determine if fasta is already in tabular format
    is_tabular = False
    standard_aa = '^[AC-IK-WY]+$'
    standard_na = '^[ACGTN]+$'
    inFile = open(fastaFile,'r')
    try:
      nseq = 0
      for i, line in enumerate( inFile ):
        line = line.rstrip( '\r\n' )
        if not line or line.startswith( '#' ):
          continue
        fields = line.split('\t')
        if len(fields) < 2: 
          is_tabular = False
          if line[0] != '>':
            print >> sys.stderr, "failed: %s does not appear to be a fasta file" % fastaFile
            exit(4)
          break
        if re.match('^[A-Z]+$',fields[-1].upper()): 
          is_tabular = True
          nseq += 1
        else:
          if line[0] != '>':
            print >> sys.stderr, "failed: %s does not appear to be a fasta file" % fastaFile
            exit(4)
        if nseq > 10:
          break
    finally:
      if inFile:
        inFile.close()
    if not is_tabular:
      fastaFile = tempfile.NamedTemporaryFile(prefix='pep_fasta_',suffix=".tab",dir=os.getcwd()).name
      fasta_to_tabular(options.protein_fasta,fastaFile)
  # vcf input 
  if options.input != None:
    try:
      inputPath = os.path.abspath(options.input)
      inputFile = open(inputPath, 'r')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(2)
  else:
    inputFile = sys.stdin
  # output 
  if options.output != None:
    try:
      outputPath = os.path.abspath(options.output)
      outputFile = open(outputPath, 'w')
    except Exception, e:
      print >> sys.stderr, "failed: %s" % e
      exit(3)
  else:
    outputFile = sys.stdout
  ## Amino_Acid_Change notations
  # G528R
  # p.Gly528Arg/c.1582G>C
  aa_change_regex = '([A-Z])(\d+)([A-Z])' # G528R
  aa_hgvs_regex = 'p\.([A-Z][a-z][a-z])(\d+)([A-Z][a-z][a-z])(/c\.(\d+)([ACGTN])>([ACGTN]))' # p.Gly528Arg/c.1582G>C
  # Save VCF file header, not currently used
  vcf_header = [] 
  reading_entries = False
  try:
    for linenum,line in enumerate(inputFile):
      ## print >> sys.stderr, "%d: %s\n" % (linenum,line)
      if line.startswith('##'):
        vcf_header.append(line)
        # May need to check SnpEff version in the header, the EFF info changed between versions 2 and 3
        ##SnpEffVersion
      elif line.startswith('#CHROM'):
        reading_entries = True
      else:
        fields = line.split('\t')
        # This is the current format of the EFF entry:
        # EFF=missense(MODERATE|MISSENSE|Ggg/Cgg|G528R|802|SCNN1D|protein_coding|CODING|ENST00000379116|12|1);OICR=(ENST00000379116|1808) 
        # If this becomes variable, will need to dynamically pattern this on the defintion in the vcf header:
        ##INFO=<ID=EFF,Number=.,Type=String,Description="Predicted effects for this variant.Format: 'Effect ( Effect_Impact | Functional_Class | Codon_Change | Amino_Acid_Change| Amino_Acid_length | Gene_Name | Transcript_BioType | Gene_Coding | Transcript_ID | Exon_Rank | Genotype_Number [ | ERRORS | WARNINGS ] )' ">
        (chrom,pos,id,ref,alts,qual,filter,info) = fields[0:8]
        for info_item in info.split(';'):
          try:
            if info_item.find('=') < 0:
              continue
            (key,val) = info_item.split('=',1)
            if key == 'EFF':
              effects = val.split(',')
              for effect in effects:
                (eff,effs) = effect.rstrip(')').split('(')
                if not (eff == 'NON_SYNONYMOUS_CODING' or eff == 'missense_variant'):
                  continue
                eff_fields = effs.split('|')
                (impact,functional_class,codon_change,aa_change,aa_len,gene_name,biotype,coding,transcript,exon) = eff_fields[0:10]
                if transcript:
                  aa_pos = None # 1-based position
                  alt_aa = '_' 
                  # parse aa_change
                  # get AA change position and alternate Animo Acid
                  sap = aa_change
                  m = re.match(aa_change_regex,aa_change)
                  if m:
                    aa_pos = int(m.groups()[1])
                    alt_aa = m.groups()[2]
                  else:
                    m = re.match(aa_hgvs_regex,aa_change)
                    if m:
                      aa_pos = int(m.groups()[1])
                      ref_aa = aa_abbrev_dict[m.groups()[0]]
                      alt_aa = aa_abbrev_dict[m.groups()[2]]
                      sap = "%s%d%s" % (ref_aa,aa_pos,alt_aa)
                  if not aa_pos:
                    continue
                  # get AA sequence
                  aa_offset = aa_pos - 1
                  (pep_id,pep_seq) = get_sequence(transcript,fastaFile)
                  if not pep_seq:
                    continue
                  start_pos = max(aa_offset - options.leading_aa_num, 0) if options.leading_aa_num else 0
                  end_pos = min(aa_offset + options.trailing_aa_num + 1, len(pep_seq)) if options.trailing_aa_num else len(pep_seq)
                  # transform sequence
                  alt_seq = pep_seq[start_pos:aa_offset] + alt_aa + pep_seq[aa_offset+1:end_pos]
                  # >ENSP00000363782 pep:known chromosome:GRCh37:1:22778472:22853855:1 gene:ENSG00000184677 transcript:ENST00000374651 gene_biotype:protein_coding transcript_biotype:protein_coding snp_location:1:22778472 codon_change:Gtg/Atg sap:V885M
                  pep_id = re.sub('pep:[a-z]*','pep:sap',pep_id)
                  hdr = ">%s snp_location:%s:%s codon_change:%s sap:%s id:%s\n" % (pep_id, chrom, pos, codon_change, sap,id)
                  outputFile.write(hdr)
                  if options.debug:
                    trimmed_seq = pep_seq[start_pos:end_pos]
                    outputFile.write(trimmed_seq)
                    outputFile.write('\n')
                  outputFile.write(alt_seq)
                  outputFile.write('\n')
          except Exception, e:
            print >> sys.stderr, "failed: %s" % e
  except Exception, e:
    print >> sys.stderr, "failed: %s" % e
    exit(1)

if __name__ == "__main__" : __main__()

