#############
# VARIABLES #
# Dangerous #
#############

# %acc2pro is used to correlate accession numbers from mascot searches to proteins
# Generally, this is added to on an as-needed basis, derived from actual mascot searches
# I have yet (6/1/2010) been able to find a consistent source for the accession numbers
# presented by mascot in the exported csv files
#
# For E coli, it's easy - multiple accession numbers can all match to the same protein
# (eg S2) as necessary.  For yeast/human, there are the S4A/S4B issues
# For the digest, there are three possible things reported:
# S4A -> peptide occurs only in the S4A variant
# S4B -> peptide occurs only in the S4B variant
# S4X -> peptide occurs in both S4A and S4B variants (most common case)
# 
# Mascot will not report S4X - it will report S4 (or maybe S4A or S4B?)
# So map all possible accession numbers to S4X in this case
# Then, in msc010203_mascotmsmsmatches.pl...
#   1) detect that the last character is "X"
#   2) search S4A then S4B then S4X for the peptide
#

%acc2pro_ecoli30S = (
# E coli 30S
#"P0AG67" => "S1",
#"P0A7V0" => "S2",
#	"A7ZHQ9" => "S2",
#"P0A7V3" => "S3",
#"P0A7V8" => "S4",
#"P0A7W1" => "S5",
#"P0A4D0" => "S6",
#"Q0TCB8" => "S7",
#"P0A7W7" => "S8",
#"P0A7X3" => "S9",
#"P0A7R5" => "S10",
#"P0A7R9" => "S11",
#"P0A7S3" => "S12",
#"P0A7S9" => "S13",
#"P0AG59" => "S14",
#"P0ADZ4" => "S15",
#"P0A7T3" => "S16",
#"P0AG63" => "S17",
#"P0A7T7" => "S18",
#"P0A7U3" => "S19",
#"P0A7U7" => "S20",
#"Q0TD41" => "S21",
"P0AG69" => "S1", # 30S ribosomal protein S1 OS=Escherichia coli O157:H7
"P0AG68" => "S1", # 30S ribosomal protein S1 OS=Escherichia coli O6
"P0AG67" => "S1", # 30S ribosomal protein S1 OS=Escherichia coli (strain K12)
"A7ZHQ9" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UIL4" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MBF0" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LGN0" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7V2" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O157:H7
"B5Z0E7" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NID0" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MP29" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O81 (strain ED1a)
"B7M1A9" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O8 (strain IAI1)
"C4ZRR1" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli (strain K12 / BW2952)
"B1XD38" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli (strain K12 / DH10B)
"A7ZWB5" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O9:H4 (strain HS)
"A1A7L3" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O1:K1 / APEC
"Q0TLG4" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7V1" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O6
"B1IQH2" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7V0" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli (strain K12)
"B7N836" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6HZE3" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli (strain SE11)
"B1LGX0" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1RG21" => "S2", # 30S ribosomal protein S2 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSK3" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK38" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS9" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4K3" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7V5" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O157:H7
"B5YTN5" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN3" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N198" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O81 (strain ED1a)
"B7M1M8" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG9" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli (strain K12 / BW2952)
"A8A5B9" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK2" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O1:K1 / APEC
"Q0TCE7" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7V4" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O6
"B1IPY5" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7V3" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli (strain K12)
"B7NDT5" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I228" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli (strain SE11)
"B1LHC8" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R612" => "S3", # 30S ribosomal protein S3 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSI5" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK20" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCR2" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LHZ8" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7W0" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O157:H7
"B5YT15" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLL6" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N181" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O81 (strain ED1a)
"B7M102" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O8 (strain IAI1)
"C4ZUF1" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli (strain K12 / BW2952)
"B1X6E8" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli (strain K12 / DH10B)
"A8A5A1" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O9:H4 (strain HS)
"A1AGI7" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O1:K1 / APEC
"Q0TCG5" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7V9" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O6
"B1IQ03" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7V8" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli (strain K12)
"B7NDR8" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I210" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli (strain SE11)
"B1LHB0" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R636" => "S4", # 30S ribosomal protein S4 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ2" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"P0A7W3" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli O157:H7
"A8A5A8" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ2" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli O1:K1 / APEC
"Q0TCF8" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7W2" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli O6
"B1IPZ6" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7W1" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli (strain K12)
"Q1R627" => "S5", # 30S ribosomal protein S5 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZV71" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UQK9" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MLK5" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LCQ6" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli (strain 55989 / EAEC)
"P0A4D1" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O157:H7
"B5Z2K6" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NTQ6" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MST0" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O81 (strain ED1a)
"B7M9G2" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O8 (strain IAI1)
"C4ZR77" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli (strain K12 / BW2952)
"B1XDV1" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli (strain K12 / DH10B)
"A8A7U6" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O9:H4 (strain HS)
"A1AJA5" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O1:K1 / APEC
"P0A4D0" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O6
"B1IT06" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P02358" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli (strain K12)
"B7NGD4" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I2A6" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli (strain SE11)
"B1LQM0" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R360" => "S6", # 30S ribosomal protein S6 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSL6" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK51" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCV6" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4L6" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli (strain 55989 / EAEC)
"P66607" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O157:H7
"B5YTP8" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLP6" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N1C2" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O81 (strain ED1a)
"B7M1P2" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O8 (strain IAI1)
"C4ZUJ6" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli (strain K12 / BW2952)
"B1X6J1" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli (strain K12 / DH10B)
"A8A5E8" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O9:H4 (strain HS)
"A1AGM8" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O1:K1 / APEC
"Q0TCB8" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P66606" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O6
"B1IPV8" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P02359" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli (strain K12)
"B7NDU9" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I241" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli (strain SE11)
"B1LHE1" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R5U2" => "S7", # 30S ribosomal protein S7 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ5" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK30" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS1" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4J4" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7W9" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O157:H7
"B5YTM7" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM5" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N190" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O81 (strain ED1a)
"B7M1M0" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG1" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli (strain K12 / BW2952)
"B1X6F8" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B1" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ5" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O1:K1 / APEC
"Q0TCF5" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7W8" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O6
"B1IPZ3" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7W7" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli (strain K12)
"B7NDS7" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I220" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli (strain SE11)
"B1LHC0" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R622" => "S8", # 30S ribosomal protein S8 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSC4" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UJW2" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MBZ1" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LHT8" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7X5" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O157:H7
"B5YSV6" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NKU2" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N102" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O81 (strain ED1a)
"B7M0U2" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O8 (strain IAI1)
"C4ZSW8" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli (strain K12 / BW2952)
"B1XHK3" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli (strain K12 / DH10B)
"A8A539" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O9:H4 (strain HS)
"A1AGC3" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O1:K1 / APEC
"Q0TCN6" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7X4" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O6
"B1IQP9" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7X3" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli (strain K12)
"B7NDK8" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I1U8" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli (strain SE11)
"B1LGJ6" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R6B0" => "S9", # 30S ribosomal protein S9 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSL0" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK45" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCT6" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4L0" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7R7" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O157:H7
"B5YTP2" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLP0" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N1A5" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O81 (strain ED1a)
"B7M1N5" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O8 (strain IAI1)
"C4ZUH6" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli (strain K12 / BW2952)
"B1X6H2" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli (strain K12 / DH10B)
"A8A5C6" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK8" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O1:K1 / APEC
"Q0TCE0" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7R6" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O6
"B1IPX8" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7R5" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli (strain K12)
"B7NDU2" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I235" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli (strain SE11)
"B1LHD5" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R601" => "S10", # 30S ribosomal protein S10 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSI6" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK21" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCR3" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LHZ9" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7S1" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O157:H7
"B5YT16" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLL7" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N182" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O81 (strain ED1a)
"B7M103" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O8 (strain IAI1)
"C4ZUF2" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli (strain K12 / BW2952)
"B1X6E9" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli (strain K12 / DH10B)
"A8A5A2" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O9:H4 (strain HS)
"A1AGI8" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O1:K1 / APEC
"Q0TCG4" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7S0" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O6
"B1IQ02" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7R9" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli (strain K12)
"B7NDR9" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I211" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli (strain SE11)
"B1LHB1" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R635" => "S11", # 30S ribosomal protein S11 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSL7" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK52" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCV7" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4L7" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7S5" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O157:H7
"B5YTP9" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLP7" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N1C3" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O81 (strain ED1a)
"B7M1P3" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O8 (strain IAI1)
"C4ZUJ7" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli (strain K12 / BW2952)
"A8A5E9" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O9:H4 (strain HS)
"A1AGM9" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O1:K1 / APEC
"Q0TCB7" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7S4" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O6
"B1IPV7" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7S3" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli (strain K12)
"B7NDV0" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I242" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli (strain SE11)
"B1LHE2" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R5U1" => "S12", # 30S ribosomal protein S12 OS=Escherichia coli (strain UTI89 / UPEC)
"P0A7T1" => "S13", # 30S ribosomal protein S13 OS=Escherichia coli O157:H7
"A1AGI9" => "S13", # 30S ribosomal protein S13 OS=Escherichia coli O1:K1 / APEC
"Q0TCG3" => "S13", # 30S ribosomal protein S13 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7T0" => "S13", # 30S ribosomal protein S13 OS=Escherichia coli O6
"P0A7S9" => "S13", # 30S ribosomal protein S13 OS=Escherichia coli (strain K12)
"Q1R633" => "S13", # 30S ribosomal protein S13 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ6" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7MCS2" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"P0AG61" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O157:H7
"B5YTM8" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM6" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7M1M1" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O8 (strain IAI1)
"B1X6F9" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B2" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ6" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O1:K1 / APEC
"Q0TCF4" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0AG60" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O6
"B1IPZ2" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0AG59" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli (strain K12)
"B7NDS8" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I221" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli (strain SE11)
"B1LHC1" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R621" => "S14", # 30S ribosomal protein S14 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZS62" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UJ60" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MB86" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LGI7" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli (strain 55989 / EAEC)
"Q8X9M2" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O157:H7
"B5YS55" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NKN4" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0V0" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O81 (strain ED1a)
"B7M073" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O8 (strain IAI1)
"C4ZSQ6" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli (strain K12 / BW2952)
"B1XGX7" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli (strain K12 / DH10B)
"A8A4Y1" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O9:H4 (strain HS)
"A1AG70" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O1:K1 / APEC
"Q0TCU4" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0ADZ5" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O6
"B1IQV6" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0ADZ4" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli (strain K12)
"B7NDF1" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I1P0" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli (strain SE11)
"B1LFR7" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R6H3" => "S15", # 30S ribosomal protein S15 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZQ49" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UH58" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MIU7" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LDJ8" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7T5" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O157:H7
"B5Z228" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NSA8" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MYP8" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O81 (strain ED1a)
"B7M978" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O8 (strain IAI1)
"C4ZYM7" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli (strain K12 / BW2952)
"B1XBT1" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli (strain K12 / DH10B)
"A8A3B6" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O9:H4 (strain HS)
"Q0TEN0" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7T4" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O6
"B1IVM4" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7T3" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli (strain K12)
"B7N6J5" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I631" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli (strain SE11)
"B1LPB6" => "S16", # 30S ribosomal protein S16 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"A7ZSK0" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK35" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS6" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4K0" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli (strain 55989 / EAEC)
"P0AG65" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O157:H7
"B5YTN2" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN0" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0V1" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O81 (strain ED1a)
"B7M1M5" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG6" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G3" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B6" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK0" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O1:K1 / APEC
"Q0TCF0" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0AG64" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O6
"B1IPY8" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0AG63" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli (strain K12)
"B7NDT2" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I225" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli (strain SE11)
"B1LHC5" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R616" => "S17", # 30S ribosomal protein S17 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZV73" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UQL1" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MLK7" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LCR3" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7T9" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O157:H7
"B5Z2K7" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NTQ7" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MST2" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O81 (strain ED1a)
"B7M9G4" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O8 (strain IAI1)
"C4ZR79" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli (strain K12 / BW2952)
"B1XDV2" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli (strain K12 / DH10B)
"A8A7U8" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O9:H4 (strain HS)
"Q0T9J1" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7T8" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O6
"B1IT04" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7T7" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli (strain K12)
"B7NGD6" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I2A8" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli (strain SE11)
"B1LQM1" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R358" => "S18", # 30S ribosomal protein S18 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSK5" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK40" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCT1" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4K5" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7U5" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O157:H7
"B5YTN7" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN5" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N1A0" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O81 (strain ED1a)
"B7M1N0" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O8 (strain IAI1)
"C4ZUH1" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G7" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli (strain K12 / DH10B)
"A8A5C1" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK4" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O1:K1 / APEC
"Q0TCE5" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7U4" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O6
"B1IPY3" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7U3" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli (strain K12)
"B7NDT7" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I230" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli (strain SE11)
"B1LHD0" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R609" => "S19", # 30S ribosomal protein S19 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZHB2" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UI68" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MAE3" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4E5" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7U9" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O157:H7
"B5YYB6" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NHC8" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MNM8" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O81 (strain ED1a)
"B7M0B9" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O8 (strain IAI1)
"C4ZPU9" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli (strain K12 / BW2952)
"B1XBE8" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli (strain K12 / DH10B)
"A7ZVX1" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O9:H4 (strain HS)
"Q0TLW7" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7U8" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O6
"B1IRF1" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7U7" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli (strain K12)
"B7N7P7" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6HZ18" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli (strain SE11)
"B1LFV4" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1RGH9" => "S20", # 30S ribosomal protein S20 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZRU7" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UIX3" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MB01" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LH01" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli (strain 55989 / EAEC)
"P68681" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O157:H7
"B5YRA5" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NJS8" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0L5" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O81 (strain ED1a)
"B7LZL5" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O8 (strain IAI1)
"C4ZQY2" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli (strain K12 / BW2952)
"B1XG70" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli (strain K12 / DH10B)
"A8A4M2" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O9:H4 (strain HS)
"Q0TD41" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P68680" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O6
"B1IRQ1" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P68679" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli (strain K12)
"B7ND54" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I437" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli (strain SE11)
"B1LF57" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R6R5" => "S21", # 30S ribosomal protein S21 OS=Escherichia coli (strain UTI89 / UPEC)
);

%acc2pro_ecoli50S = (
# E coli 50S
#"P0A7L0" => "L1",
#"Q0TCE4" => "L2",
#"Q0TCE2" => "L4",
#"Q0TCF3" => "L5",
#"P0AG55" => "L6",
#"P0A7K2" => "L7L12",
#"P0A7R1" => "L9",
#"P0A7J3" => "L10",
#"P0A7J7" => "L11",
#"Q0TCG0" => "L15",
#"P0ADY7" => "L16",
#"P0AG44" => "L17",
#"P0C018" => "L18",
#"P0AG48" => "L21",
#"B1IY84" => "L25",
#"P0A7L8" => "L27",
#"P0A7M2" => "L28",
#"P0A7M6" => "L29",
#"P0AG51" => "L30",
#"P0A7M9" => "L31",
#"P0A7N9" => "L33",
"A7ZUJ7" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UPD9" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MIX0" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LA77" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7L2" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O157:H7
"B5Z080" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NRR2" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MR70" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O81 (strain ED1a)
"B7M731" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O8 (strain IAI1)
"C5A0S4" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli (strain K12 / BW2952)
"B1XBY6" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli (strain K12 / DH10B)
"A8A783" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O9:H4 (strain HS)
"Q0TA81" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7L1" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O6
"B1IUR3" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7L0" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli (strain K12)
"B7NFS4" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I5J4" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli (strain SE11)
"B1LNT6" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R5U9" => "L1", # 50S ribosomal protein L1 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSK6" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK41" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCT2" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4K6" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli (strain 55989 / EAEC)
"P60424" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O157:H7
"B5YTN8" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN6" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N1A1" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O81 (strain ED1a)
"B7M1N1" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O8 (strain IAI1)
"C4ZUH2" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G8" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli (strain K12 / DH10B)
"A8A5C2" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK5" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O1:K1 / APEC
"Q0TCE4" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P60423" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O6
"B1IPY2" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P60422" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli (strain K12)
"B7NDT8" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I231" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli (strain SE11)
"B1LHD1" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R607" => "L2", # 50S ribosomal protein L2 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSK9" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK44" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCT5" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4K9" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli (strain 55989 / EAEC)
"P60440" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O157:H7
"B5YTP1" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN9" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0W0" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O81 (strain ED1a)
"B7M1N4" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O8 (strain IAI1)
"C4ZUH5" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli (strain K12 / BW2952)
"B1X6H1" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli (strain K12 / DH10B)
"A8A5C5" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK7" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O1:K1 / APEC
"Q0TCE1" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P60439" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O6
"B1IPX9" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P60438" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli (strain K12)
"B7NDU1" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I234" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli (strain SE11)
"B1LHD4" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R602" => "L3", # 50S ribosomal protein L3 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSK8" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK43" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCT4" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4K8" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli (strain 55989 / EAEC)
"P60725" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O157:H7
"B5YTP0" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN8" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N1A3" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O81 (strain ED1a)
"B7M1N3" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O8 (strain IAI1)
"C4ZUH4" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli (strain K12 / BW2952)
"B1X6H0" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli (strain K12 / DH10B)
"A8A5C4" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK6" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O1:K1 / APEC
"Q0TCE2" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P60724" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O6
"B1IPY0" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P60723" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli (strain K12)
"B7NDU0" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I233" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli (strain SE11)
"B1LHD3" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R604" => "L4", # 50S ribosomal protein L4 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ7" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK32" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS3" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4J7" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli (strain 55989 / EAEC)
"P62401" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O157:H7
"B5YTM9" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM7" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N192" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O81 (strain ED1a)
"B7M1M2" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG3" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G0" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B3" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ7" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O1:K1 / APEC
"Q0TCF3" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P62400" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O6
"B1IPZ1" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P62399" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli (strain K12)
"B7NDS9" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I222" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli (strain SE11)
"B1LHC2" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R620" => "L5", # 50S ribosomal protein L5 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ4" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK29" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS0" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LI04" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli (strain 55989 / EAEC)
"P0AG57" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O157:H7
"B5YTM6" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM4" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N189" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O81 (strain ED1a)
"B7M1L9" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG0" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli (strain K12 / BW2952)
"B1X6F7" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B0" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ4" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O1:K1 / APEC
"Q0TCF6" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0AG56" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O6
"B1IPZ4" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0AG55" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli (strain K12)
"B7NDS6" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I219" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli (strain SE11)
"B1LHB9" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R624" => "L6", # 50S ribosomal protein L6 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZUK0" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UPE1" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MIX2" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LA79" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7K4" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O157:H7
"B5Z082" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NRR4" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MRB2" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O81 (strain ED1a)
"B7M733" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O8 (strain IAI1)
"C5A0S6" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli (strain K12 / BW2952)
"B1XBY8" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli (strain K12 / DH10B)
"A8A785" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O9:H4 (strain HS)
"Q0TA79" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7K3" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O6
"B1IUR1" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7K2" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli (strain K12)
"B7NFS6" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I5J6" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli (strain SE11)
"B1LNT8" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R5V1" => "L7L12", # 50S ribosomal protein L7/L12 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZV74" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UQL2" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MLK8" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LCR4" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7R3" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O157:H7
"B5Z2K8" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NTQ8" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MT77" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O81 (strain ED1a)
"B7M9G5" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O8 (strain IAI1)
"C4ZR80" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli (strain K12 / BW2952)
"B1XDV3" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli (strain K12 / DH10B)
"A8A7U9" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O9:H4 (strain HS)
"A1AJA7" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O1:K1 / APEC
"Q0T9J0" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7R2" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O6
"B1IT03" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7R1" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli (strain K12)
"B7NGD7" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I2A9" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli (strain SE11)
"B1LR79" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R357" => "L9", # 50S ribosomal protein L9 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZUJ8" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UPE0" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MIX1" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LA78" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7J5" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O157:H7
"B5Z081" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NRR3" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MR71" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O81 (strain ED1a)
"B7M732" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O8 (strain IAI1)
"C5A0S5" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli (strain K12 / BW2952)
"B1XBY7" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli (strain K12 / DH10B)
"A8A784" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O9:H4 (strain HS)
"A1AIF7" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O1:K1 / APEC
"Q0TA80" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7J4" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O6
"B1IUR2" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7J3" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli (strain K12)
"B7NFS5" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I5J5" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli (strain SE11)
"B1LNT7" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R5V0" => "L10", # 50S ribosomal protein L10 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZUJ6" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UPD8" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MIW9" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LA74" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7J9" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O157:H7
"B5Z079" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NRR1" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MRA9" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O81 (strain ED1a)
"B7M730" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O8 (strain IAI1)
"C5A0S3" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli (strain K12 / BW2952)
"B1XBY5" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli (strain K12 / DH10B)
"A8A782" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O9:H4 (strain HS)
"Q0TA82" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7J8" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O6
"B1IUR4" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7J7" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli (strain K12)
"B7NFS3" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I5J3" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli (strain SE11)
"B1LNT5" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R5U7" => "L11", # 50S ribosomal protein L11 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSC5" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UJW3" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MBZ2" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LHT9" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli (strain 55989 / EAEC)
"P0AA12" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O157:H7
"B5YSV7" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NKU3" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0L6" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O81 (strain ED1a)
"B7M0U3" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O8 (strain IAI1)
"C4ZSW9" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli (strain K12 / BW2952)
"B1XHK4" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli (strain K12 / DH10B)
"A8A540" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O9:H4 (strain HS)
"A1AGC4" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O1:K1 / APEC
"Q0TCN5" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0AA11" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O6
"B1IQP8" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0AA10" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli (strain K12)
"B7NDK9" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I1U9" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli (strain SE11)
"B1LGJ7" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R6A9" => "L13", # 50S ribosomal protein L13 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ9" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK34" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS5" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4J9" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli (strain 55989 / EAEC)
"P0ADY5" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O157:H7
"B5YTN1" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM9" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N194" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O81 (strain ED1a)
"B7M1M4" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG5" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G2" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B5" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ9" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O1:K1 / APEC
"Q0TCF1" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0ADY4" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O6
"B1IPY9" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0ADY3" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli (strain K12)
"B7NDT1" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I224" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli (strain SE11)
"B1LHC4" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R617" => "L14", # 50S ribosomal protein L14 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ0" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK25" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCR6" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LI02" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli (strain 55989 / EAEC)
"P66072" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O157:H7
"B5YTM2" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM0" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N185" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O81 (strain ED1a)
"B7M106" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O8 (strain IAI1)
"C4ZUF6" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli (strain K12 / BW2952)
"B1X6F3" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli (strain K12 / DH10B)
"A8A5A6" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ1" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O1:K1 / APEC
"Q0TCG0" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P66071" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O6
"B1IPZ8" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P02413" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli (strain K12)
"B7NDS2" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I215" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli (strain SE11)
"B1LHB5" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R630" => "L15", # 50S ribosomal protein L15 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSK2" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK37" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS8" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4K2" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli (strain 55989 / EAEC)
"P0ADY9" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O157:H7
"B5YTN4" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN2" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N197" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O81 (strain ED1a)
"B7M1M7" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG8" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G5" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B8" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK1" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O1:K1 / APEC
"Q0TCE8" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0ADY8" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O6
"B1IPY6" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0ADY7" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli (strain K12)
"B7NDT4" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I227" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli (strain SE11)
"B1LHC7" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R613" => "L16", # 50S ribosomal protein L16 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSI3" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK18" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCR0" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LHZ6" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli (strain 55989 / EAEC)
"P0AG46" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O157:H7
"B5YT13" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLL4" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N179" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O81 (strain ED1a)
"B7M100" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O8 (strain IAI1)
"C4ZUE9" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli (strain K12 / BW2952)
"B1X6E6" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli (strain K12 / DH10B)
"A8A599" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O9:H4 (strain HS)
"A1AGI5" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O1:K1 / APEC
"Q0TCG7" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0AG45" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O6
"B1IQ05" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0AG44" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli (strain K12)
"B7NDR6" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I208" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli (strain SE11)
"B1LGQ0" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R638" => "L17", # 50S ribosomal protein L17 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ3" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK28" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCR9" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LI05" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli (strain 55989 / EAEC)
"P0C020" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O157:H7
"B5YTM5" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM3" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0U4" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O81 (strain ED1a)
"B7M1L8" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O8 (strain IAI1)
"C4ZUF9" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli (strain K12 / BW2952)
"B1X6F6" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli (strain K12 / DH10B)
"A8A5A9" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ3" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O1:K1 / APEC
"Q0TCF7" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0C019" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O6
"B1IPZ5" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0C018" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli (strain K12)
"B7NDS5" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I218" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli (strain SE11)
"B1LHB8" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R626" => "L18", # 50S ribosomal protein L18 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZQ46" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UH55" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MIU4" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LDJ5" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7K8" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O157:H7
"B5Z225" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NSA5" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MYP6" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O81 (strain ED1a)
"B7M975" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O8 (strain IAI1)
"C4ZYM4" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli (strain K12 / BW2952)
"B1XBS8" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli (strain K12 / DH10B)
"A8A3B3" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O9:H4 (strain HS)
"A1AED6" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O1:K1 / APEC
"Q0TEN3" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7K7" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O6
"B1IVM7" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7K6" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli (strain K12)
"B7N6J2" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I628" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli (strain SE11)
"B1LPB3" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R8B9" => "L19", # 50S ribosomal protein L19 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZMI4" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7US99" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MAS6" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L6I9" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7L5" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O157:H7
"B5YQ04" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NT61" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MVJ4" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O81 (strain ED1a)
"B7M1C4" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O8 (strain IAI1)
"C4ZYH8" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli (strain K12 / BW2952)
"B1XG23" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli (strain K12 / DH10B)
"A8A0Q9" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O9:H4 (strain HS)
"A1ABQ0" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O1:K1 / APEC
"Q0THB3" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7L4" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O6
"B1IPL2" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7L3" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli (strain K12)
"B7N553" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6IBD4" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli (strain SE11)
"B1LE15" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1RB79" => "L20", # 50S ribosomal protein L20 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZS83" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UJ79" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MBV5" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LHP9" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli (strain 55989 / EAEC)
"P0AG49" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O157:H7
"B5YS75" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NKQ3" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0W8" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O81 (strain ED1a)
"B7M092" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O8 (strain IAI1)
"C4ZSS5" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli (strain K12 / BW2952)
"B1XHG2" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli (strain K12 / DH10B)
"A8A501" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O9:H4 (strain HS)
"Q0TCS4" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"B1IQT7" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0AG48" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli (strain K12)
"B7NDH0" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I1Q9" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli (strain SE11)
"B1LGF2" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R6F0" => "L21", # 50S ribosomal protein L21 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSK4" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK39" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCT0" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4K4" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli (strain 55989 / EAEC)
"P61177" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O157:H7
"B5YTN6" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN4" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N199" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O81 (strain ED1a)
"B7M1M9" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O8 (strain IAI1)
"C4ZUH0" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G6" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli (strain K12 / DH10B)
"A8A5C0" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O9:H4 (strain HS)
"A1AGK3" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O1:K1 / APEC
"Q0TCE6" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P61176" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O6
"B1IPY4" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P61175" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli (strain K12)
"B7NDT6" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I229" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli (strain SE11)
"B1LHC9" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R610" => "L22", # 50S ribosomal protein L22 OS=Escherichia coli (strain UTI89 / UPEC)
"P0ADZ2" => "L23", # 50S ribosomal protein L23 OS=Escherichia coli O157:H7
"Q0TCE3" => "L23", # 50S ribosomal protein L23 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0ADZ1" => "L23", # 50S ribosomal protein L23 OS=Escherichia coli O6
"P0ADZ0" => "L23", # 50S ribosomal protein L23 OS=Escherichia coli (strain K12)
"Q1R606" => "L23", # 50S ribosomal protein L23 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ8" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK33" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS4" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4J5" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli (strain 55989 / EAEC)
"P60625" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O157:H7
"B5YTN0" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM8" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N193" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O81 (strain ED1a)
"B7M1M3" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG4" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G1" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B4" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O9:H4 (strain HS)
"A1AGJ8" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O1:K1 / APEC
"Q0TCF2" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"Q8FD03" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O6
"B1IPZ0" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P60624" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli (strain K12)
"B7NDT0" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I223" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli (strain SE11)
"B1LHC3" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R619" => "L24", # 50S ribosomal protein L24 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZP09" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UFK2" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MFA0" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LAK9" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli (strain 55989 / EAEC)
"Q8XE69" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O157:H7
"B5YWX8" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NN00" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MXK1" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O81 (strain ED1a)
"B7M535" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O8 (strain IAI1)
"C4ZU30" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli (strain K12 / BW2952)
"B1X883" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli (strain K12 / DH10B)
"A8A248" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O9:H4 (strain HS)
"Q0TFQ4" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"Q8FFS1" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O6
"B1IY84" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P68919" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli (strain K12)
"B7N5E9" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I184" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli (strain SE11)
"B1LKT4" => "L25", # 50S ribosomal protein L25 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"A7ZS81" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UJ78" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MBV4" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LHP5" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7M0" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O157:H7
"B5YS74" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NKQ2" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0H8" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O81 (strain ED1a)
"B7M091" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O8 (strain IAI1)
"C4ZSS4" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli (strain K12 / BW2952)
"B1XHG1" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli (strain K12 / DH10B)
"A8A500" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O9:H4 (strain HS)
"Q0TCS5" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7L9" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O6
"B1IQT8" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7L8" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli (strain K12)
"B7NDG9" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I1Q8" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli (strain SE11)
"B1LGF1" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R6F2" => "L27", # 50S ribosomal protein L27 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZTI8" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7ULJ2" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MFJ8" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L760" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7M4" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O157:H7
"B5YWD5" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NPZ8" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N277" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O81 (strain ED1a)
"B7M4C1" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O8 (strain IAI1)
"C4ZXM9" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli (strain K12 / BW2952)
"B1X971" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli (strain K12 / DH10B)
"A8A699" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O9:H4 (strain HS)
"Q0TBH2" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7M3" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O6
"B1IZF6" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7M2" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli (strain K12)
"B7NEU1" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I3L4" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli (strain SE11)
"B1LK74" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R4V5" => "L28", # 50S ribosomal protein L28 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSK1" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK36" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCS7" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L4J8" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7M8" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O157:H7
"B5YTN3" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLN1" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N196" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O81 (strain ED1a)
"B7M1M6" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O8 (strain IAI1)
"C4ZUG7" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli (strain K12 / BW2952)
"B1X6G4" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli (strain K12 / DH10B)
"A8A5B7" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O9:H4 (strain HS)
"Q0TCE9" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7M7" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O6
"B1IPY7" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7M6" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli (strain K12)
"B7NDT3" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I226" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli (strain SE11)
"B1LHC6" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R615" => "L29", # 50S ribosomal protein L29 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSJ1" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UK26" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCR7" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LI03" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli (strain 55989 / EAEC)
"P0AG53" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O157:H7
"B5YTM3" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NLM1" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N0U2" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O81 (strain ED1a)
"B7M107" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O8 (strain IAI1)
"C4ZUF7" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli (strain K12 / BW2952)
"B1X6F4" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli (strain K12 / DH10B)
"A8A5A7" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O9:H4 (strain HS)
"Q0TCF9" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0AG52" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O6
"B1IPZ7" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0AG51" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli (strain K12)
"B7NDS3" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I216" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli (strain SE11)
"B1LHB6" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R629" => "L30", # 50S ribosomal protein L30 OS=Escherichia coli (strain UTI89 / UPEC)
"P0A7N3" => "L31", # 50S ribosomal protein L31 type B 1 OS=Escherichia coli O157:H7
"Q8X9T8" => "L31", # 50S ribosomal protein L31 type B 2 OS=Escherichia coli O157:H7
"B7UJE5" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MCA6" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L427" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli (strain 55989 / EAEC)
"B7NK70" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7M2U3" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli O8 (strain IAI1)
"B1XE39" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli (strain K12 / DH10B)
"A7ZWT4" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli O9:H4 (strain HS)
"P0A7N2" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli O6
"P0A7N1" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli (strain K12)
"B7N8J3" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I084" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli (strain SE11)
"B1LHW4" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1RFP9" => "L31", # 50S ribosomal protein L31 type B OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZUF1" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UNQ6" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MI68" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LA33" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7N0" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O157:H7
"B5YZ77" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NU69" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N2S7" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O81 (strain ED1a)
"B7M6Y7" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O8 (strain IAI1)
"C5A0A3" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli (strain K12 / BW2952)
"B1XBA2" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli (strain K12 / DH10B)
"A8A742" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O9:H4 (strain HS)
"Q0TAC8" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0C076" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O6
"B1IVE3" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7M9" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli (strain K12)
"B7NFN5" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I4S9" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli (strain SE11)
"B1LNP1" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"P0C203" => "L31", # 50S ribosomal protein L31 OS=Escherichia coli (strain UTI89 / UPEC)
"B7UPA5" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MJ76" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7LG22" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7N5" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O157:H7
"B5YVV7" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NKJ1" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MTL9" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O81 (strain ED1a)
"B7LX23" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O8 (strain IAI1)
"C4ZS29" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli (strain K12 / BW2952)
"B1X9Z9" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli (strain K12 / DH10B)
"A7ZZ46" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O9:H4 (strain HS)
"Q0TIY5" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"B1IV14" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7N4" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli (strain K12)
"B7NAW5" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I9G6" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli (strain SE11)
"B1LI56" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1RD68" => "L32", # 50S ribosomal protein L32 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZTI7" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7ULJ1" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MFJ7" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L759" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7P1" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O157:H7
"B5YWD4" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NPE2" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N276" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O81 (strain ED1a)
"B7M4C0" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O8 (strain IAI1)
"C4ZXM8" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli (strain K12 / BW2952)
"B1X970" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli (strain K12 / DH10B)
"A8A698" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O9:H4 (strain HS)
"Q0TBH3" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"P0A7P0" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O6
"B1IZF7" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7N9" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli (strain K12)
"B7NEU0" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I3L3" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli (strain SE11)
"B1LK73" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R4V6" => "L33", # 50S ribosomal protein L33 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZTQ9" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B7UMH0" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MGC4" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L847" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7P7" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O157:H7
"B5YXA7" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NR05" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7N2E8" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O81 (strain ED1a)
"B7M555" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O8 (strain IAI1)
"C4ZYY0" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli (strain K12 / BW2952)
"B1X9T1" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli (strain K12 / DH10B)
"A8A6G4" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O9:H4 (strain HS)
"P0A7P6" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O6
"B1IX36" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7P5" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli (strain K12)
"B7NF20" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6I3T6" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli (strain SE11)
"B1LL30" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R4N3" => "L34", # 50S ribosomal protein L34 OS=Escherichia coli (strain UTI89 / UPEC)
"B7USA0" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O127:H6 (strain E2348/69 / EPEC)
"B7MAS7" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O45:K1 (strain S88 / ExPEC)
"B7L6J0" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli (strain 55989 / EAEC)
"P0A7Q2" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O157:H7
"B5YQ05" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O157:H7 (strain EC4115 / EHEC)
"B7NT60" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"B7MVJ5" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O81 (strain ED1a)
"B7M1C5" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O8 (strain IAI1)
"C4ZYH9" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli (strain K12 / BW2952)
"B1XG24" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli (strain K12 / DH10B)
"A8A0R0" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O9:H4 (strain HS)
"Q0THB2" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"B1IPL1" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"P0A7Q1" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli (strain K12)
"B7N554" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli O17:K52:H18 (strain UMN026 / ExPEC)
"B6IBD5" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli (strain SE11)
"B1LE14" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1RB78" => "L35", # 50S ribosomal protein L35 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZSI8" => "L36", # 50S ribosomal protein L36 1 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B1X6F1" => "L36", # 50S ribosomal protein L36 1 OS=Escherichia coli (strain K12 / DH10B)
"A8A5A4" => "L36", # 50S ribosomal protein L36 1 OS=Escherichia coli O9:H4 (strain HS)
"Q0TCG2" => "L36", # 50S ribosomal protein L36 1 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"B1IQ00" => "L36", # 50S ribosomal protein L36 1 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"B1LHB3" => "L36", # 50S ribosomal protein L36 1 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1R632" => "L36", # 50S ribosomal protein L36 1 OS=Escherichia coli (strain UTI89 / UPEC)
"A7ZI30" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli O139:H28 (strain E24377A / ETEC)
"B1XE38" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli (strain K12 / DH10B)
"A7ZWT3" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli O9:H4 (strain HS)
"Q0TKY8" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli O6:K15:H31 (strain 536 / UPEC)
"Q8FKL1" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli O6
"B1J0X7" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli (strain ATCC 8739 / DSM 1576 / Crooks)
"Q2EEQ2" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli (strain K12)
"B1LHW3" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli (strain SMS-3-5 / SECEC)
"Q1RFQ0" => "L36", # 50S ribosomal protein L36 2 OS=Escherichia coli (strain UTI89 / UPEC)
"P0A7Q7" => "L36", # 50S ribosomal protein L36 OS=Escherichia coli O157:H7
"B7NK71" => "L36", # 50S ribosomal protein L36 OS=Escherichia coli O7:K1 (strain IAI39 / ExPEC)
"C4ZUF4" => "L36", # 50S ribosomal protein L36 OS=Escherichia coli (strain K12 / BW2952)
"P0A7Q6" => "L36", # 50S ribosomal protein L36 OS=Escherichia coli (strain K12)
);

%acc2pro_yeast40S = (
# Yeast 40S
# S2, S3, S5, S12, S13, S15, S20, S31
#"P32905" => "S0X",
#"P33442" => "S1X",
#"P23248" => "S1X",
#"P25443" => "S2",
#"P05750" => "S3",
#"P05753" => "S4X",
#"P26783" => "S5",
#"P02365" => "S6X",
#"P48164" => "S7X",
#"P05754" => "S8X",
#"O13516" => "S9X",
#"P05755" => "S9X",
#"Q08745" => "S10X",
#"P26781" => "S11X",
#"P48589" => "S12",
#"P05756" => "S13",
#"P06367" => "S14X",
#"Q01855" => "S15",
#"P40213" => "S16X",
#"P02407" => "S17X",
#"P35271" => "S18X",
#"P07280" => "S19X",
#"P38701" => "S20",
#"P0C0V8" => "S21X",
#"P0C0W1" => "S22X",
#"P26782" => "S24X",
#"P0C0T4" => "S25X",
#"P39938" => "S26X",
#"P35997" => "S27X",
#"P0C0X0" => "S28X",
#"P41058" => "S29X",
#"Q12087" => "S30X",
#"P05759" => "S31",

"Q08745" => "S10X", # 40S ribosomal protein S10-A OS=Saccharomyces cerevisiae
"P46784" => "S10X", # 40S ribosomal protein S10-B OS=Saccharomyces cerevisiae
"P26781" => "S11", # 40S ribosomal protein S11 OS=Saccharomyces cerevisiae
"P48589" => "S12", # 40S ribosomal protein S12 OS=Saccharomyces cerevisiae
"P05756" => "S13", # 40S ribosomal protein S13 OS=Saccharomyces cerevisiae
"P06367" => "S14X", # 40S ribosomal protein S14-A OS=Saccharomyces cerevisiae
"P39516" => "S14X", # 40S ribosomal protein S14-B OS=Saccharomyces cerevisiae
"Q01855" => "S15", # 40S ribosomal protein S15 OS=Saccharomyces cerevisiae
"P40213" => "S16", # 40S ribosomal protein S16 OS=Saccharomyces cerevisiae
"P02407" => "S17X", # 40S ribosomal protein S17-A OS=Saccharomyces cerevisiae
"P14127" => "S17X", # 40S ribosomal protein S17-B OS=Saccharomyces cerevisiae
"P35271" => "S18", # 40S ribosomal protein S18 OS=Saccharomyces cerevisiae
"P07280" => "S19X", # 40S ribosomal protein S19-A OS=Saccharomyces cerevisiae
"P07281" => "S19X", # 40S ribosomal protein S19-B OS=Saccharomyces cerevisiae
"P38701" => "S20", # 40S ribosomal protein S20 OS=Saccharomyces cerevisiae
"P0C0V8" => "S21X", # 40S ribosomal protein S21-A OS=Saccharomyces cerevisiae
"Q3E754" => "S21X", # 40S ribosomal protein S21-B OS=Saccharomyces cerevisiae
"P0C0W1" => "S22X", # 40S ribosomal protein S22-A OS=Saccharomyces cerevisiae
"Q3E7Y3" => "S22X", # 40S ribosomal protein S22-B OS=Saccharomyces cerevisiae
"P32827" => "S23", # 40S ribosomal protein S23 OS=Saccharomyces cerevisiae
"P26782" => "S24", # 40S ribosomal protein S24 OS=Saccharomyces cerevisiae
"Q3E792" => "S25X", # 40S ribosomal protein S25-A OS=Saccharomyces cerevisiae
"P0C0T4" => "S25X", # 40S ribosomal protein S25-B OS=Saccharomyces cerevisiae
"P39938" => "S26X", # 40S ribosomal protein S26-A OS=Saccharomyces cerevisiae
"P39939" => "S26X", # 40S ribosomal protein S26-B OS=Saccharomyces cerevisiae
"P35997" => "S27X", # 40S ribosomal protein S27-A OS=Saccharomyces cerevisiae
"P38711" => "S27X", # 40S ribosomal protein S27-B OS=Saccharomyces cerevisiae
"Q3E7X9" => "S28X", # 40S ribosomal protein S28-A OS=Saccharomyces cerevisiae
"P0C0X0" => "S28X", # 40S ribosomal protein S28-B OS=Saccharomyces cerevisiae
"P41057" => "S29X", # 40S ribosomal protein S29-A OS=Saccharomyces cerevisiae
"P41058" => "S29X", # 40S ribosomal protein S29-B OS=Saccharomyces cerevisiae
"P25443" => "S2", # 40S ribosomal protein S2 OS=Saccharomyces cerevisiae
"Q12087" => "S30", # 40S ribosomal protein S30 OS=Saccharomyces cerevisiae
"P05759" => "S31", # 40S ribosomal protein S31 OS=Saccharomyces cerevisiae
"B3RHV0" => "S1X", # 40S ribosomal protein S1-A OS=Saccharomyces cerevisiae (strain RM11-1a)
"C7GSL4" => "S1X", # 40S ribosomal protein S1-A OS=Saccharomyces cerevisiae (strain JAY291)
"B5VNY0" => "S1X", # 40S ribosomal protein S1-A OS=Saccharomyces cerevisiae (strain AWRI1631)
"A7A1W1" => "S1X", # 40S ribosomal protein S1-A OS=Saccharomyces cerevisiae (strain YJM789)
"P33442" => "S1X", # 40S ribosomal protein S1-A OS=Saccharomyces cerevisiae
"B3LLJ2" => "S1X", # 40S ribosomal protein S1-B OS=Saccharomyces cerevisiae (strain RM11-1a)
"C7GVZ9" => "S1X", # 40S ribosomal protein S1-B OS=Saccharomyces cerevisiae (strain JAY291)
"B5VP69" => "S1X", # 40S ribosomal protein S1-B OS=Saccharomyces cerevisiae (strain AWRI1631)
"A6ZM02" => "S1X", # 40S ribosomal protein S1-B OS=Saccharomyces cerevisiae (strain YJM789)
"P23248" => "S1X", # 40S ribosomal protein S1-B OS=Saccharomyces cerevisiae
"P05750" => "S3", # 40S ribosomal protein S3 OS=Saccharomyces cerevisiae
"P05753" => "S4", # 40S ribosomal protein S4 OS=Saccharomyces cerevisiae
"P26783" => "S5", # 40S ribosomal protein S5 OS=Saccharomyces cerevisiae
"P02365" => "S6", # 40S ribosomal protein S6 OS=Saccharomyces cerevisiae
"P26786" => "S7X", # 40S ribosomal protein S7-A OS=Saccharomyces cerevisiae
"P48164" => "S7X", # 40S ribosomal protein S7-B OS=Saccharomyces cerevisiae
"P05754" => "S8", # 40S ribosomal protein S8 OS=Saccharomyces cerevisiae
"O13516" => "S9X", # 40S ribosomal protein S9-A OS=Saccharomyces cerevisiae
"P05755" => "S9X", # 40S ribosomal protein S9-B OS=Saccharomyces cerevisiae
"B3LI22" => "S0X", # 40S ribosomal protein S0-A OS=Saccharomyces cerevisiae (strain RM11-1a)
"A6ZUM5" => "S0X", # 40S ribosomal protein S0-A OS=Saccharomyces cerevisiae (strain YJM789)
"P32905" => "S0X", # 40S ribosomal protein S0-A OS=Saccharomyces cerevisiae
"B3LT19" => "S0X", # 40S ribosomal protein S0-B OS=Saccharomyces cerevisiae (strain RM11-1a)
"A7A0V3" => "S0X", # 40S ribosomal protein S0-B OS=Saccharomyces cerevisiae (strain YJM789)
"P46654" => "S0X", # 40S ribosomal protein S0-B OS=Saccharomyces cerevisiae

);

%acc2pro_yeast60S = (
# Yeast 60S
# These are NOT A/B/X
# L3, L5, L10, L25, L28, L29, L30, L32, L38, L39, P0
#"P53030" => "L1X",
#"P05736" => "L2X",
#"P14126" => "L3",
#"P10664" => "L4X",
#"P26321" => "L5",
#"P05739" => "L6X",
#	"Q02326" => "L6X",
#"P05737" => "L7X",
#	"Q12213" => "L7X",
#"P17076" => "L8X",
#"P05738" => "L9X",
#	"P51401" => "L9X",
#"P41805" => "L10",
#"P0C0W9" => "L11X",
#"Q12690" => "L13X",
#"P36105" => "L14X",
#"P05748" => "L15X",
#"P26785" => "L16X",
#	"P26784" => "L16X",
#"P05740" => "L17X",
#	"P46990" => "L17X",
#"P07279" => "L18X",
#"P05735" => "L19X",
#"P0C2I0" => "L20X",
#"Q02753" => "L21X",
#"P04451" => "L23X",
#"P04449" => "L24X",
#"P24000" => "L24X",
#"P04456" => "L25",
#"P05743" => "L26X",
#	"P53221" => "L26X",
#"P0C2H6" => "L27X",
#"P02406" => "L28",
#"P14120" => "L30",
#"P0C2H8" => "L31X",
#"P38061" => "L32",
#"P05744" => "L33X",
#"P41056" => "L33X",
#"P39741" => "L35X",
#"O14455" => "L36X",
#	"P05745" => "L36X",
#"P49166" => "L37X",
#"P49167" => "L38",
#"P02405" => "L42X",
#"P49631" => "L43X",
#"P05317" => "P0",
#"P02400" => "P2",

"P41805" => "L10", # 60S ribosomal protein L10 OS=Saccharomyces cerevisiae
"P0C0W9" => "L11X", # 60S ribosomal protein L11-A OS=Saccharomyces cerevisiae
"Q3E757" => "L11X", # 60S ribosomal protein L11-B OS=Saccharomyces cerevisiae
"P17079" => "L12", # 60S ribosomal protein L12 OS=Saccharomyces cerevisiae
"Q12690" => "L13X", # 60S ribosomal protein L13-A OS=Saccharomyces cerevisiae
"P40212" => "L13X", # 60S ribosomal protein L13-B OS=Saccharomyces cerevisiae
"P36105" => "L14X", # 60S ribosomal protein L14-A OS=Saccharomyces cerevisiae
"P38754" => "L14X", # 60S ribosomal protein L14-B OS=Saccharomyces cerevisiae
"P05748" => "L15X", # 60S ribosomal protein L15-A OS=Saccharomyces cerevisiae
"P54780" => "L15X", # 60S ribosomal protein L15-B OS=Saccharomyces cerevisiae
"P26784" => "L16X", # 60S ribosomal protein L16-A OS=Saccharomyces cerevisiae
"P26785" => "L16X", # 60S ribosomal protein L16-B OS=Saccharomyces cerevisiae
"P05740" => "L17X", # 60S ribosomal protein L17-A OS=Saccharomyces cerevisiae
"P46990" => "L17X", # 60S ribosomal protein L17-B OS=Saccharomyces cerevisiae
"P07279" => "L18", # 60S ribosomal protein L18 OS=Saccharomyces cerevisiae
"P05735" => "L19", # 60S ribosomal protein L19 OS=Saccharomyces cerevisiae
"P53030" => "L1", # 60S ribosomal protein L1 OS=Saccharomyces cerevisiae
"P0C2I0" => "L20", # 60S ribosomal protein L20 OS=Saccharomyces cerevisiae
"Q02753" => "L21X", # 60S ribosomal protein L21-A OS=Saccharomyces cerevisiae
"Q12672" => "L21X", # 60S ribosomal protein L21-B OS=Saccharomyces cerevisiae
"P05749" => "L22X", # 60S ribosomal protein L22-A OS=Saccharomyces cerevisiae
"P56628" => "L22X", # 60S ribosomal protein L22-B OS=Saccharomyces cerevisiae
"P04451" => "L23", # 60S ribosomal protein L23 OS=Saccharomyces cerevisiae
"P04449" => "L24X", # 60S ribosomal protein L24-A OS=Saccharomyces cerevisiae
"P24000" => "L24X", # 60S ribosomal protein L24-B OS=Saccharomyces cerevisiae
"P04456" => "L25", # 60S ribosomal protein L25 OS=Saccharomyces cerevisiae
"P05743" => "L26X", # 60S ribosomal protein L26-A OS=Saccharomyces cerevisiae
"P53221" => "L26X", # 60S ribosomal protein L26-B OS=Saccharomyces cerevisiae
"P0C2H6" => "L27X", # 60S ribosomal protein L27-A OS=Saccharomyces cerevisiae
"P0C2H7" => "L27X", # 60S ribosomal protein L27-B OS=Saccharomyces cerevisiae
"P02406" => "L28", # 60S ribosomal protein L28 OS=Saccharomyces cerevisiae
"P05747" => "L29", # 60S ribosomal protein L29 OS=Saccharomyces cerevisiae
"P05736" => "L2", # 60S ribosomal protein L2 OS=Saccharomyces cerevisiae
"P14120" => "L30", # 60S ribosomal protein L30 OS=Saccharomyces cerevisiae
"P0C2H8" => "L31X", # 60S ribosomal protein L31-A OS=Saccharomyces cerevisiae
"P0C2H9" => "L31X", # 60S ribosomal protein L31-B OS=Saccharomyces cerevisiae
"P38061" => "L32", # 60S ribosomal protein L32 OS=Saccharomyces cerevisiae
"P05744" => "L33X", # 60S ribosomal protein L33-A OS=Saccharomyces cerevisiae
"P41056" => "L33X", # 60S ribosomal protein L33-B OS=Saccharomyces cerevisiae
"P87262" => "L34X", # 60S ribosomal protein L34-A OS=Saccharomyces cerevisiae
"P40525" => "L34X", # 60S ribosomal protein L34-B OS=Saccharomyces cerevisiae
"P39741" => "L35", # 60S ribosomal protein L35 OS=Saccharomyces cerevisiae
"P05745" => "L36X", # 60S ribosomal protein L36-A OS=Saccharomyces cerevisiae
"O14455" => "L36X", # 60S ribosomal protein L36-B OS=Saccharomyces cerevisiae
"P49166" => "L37X", # 60S ribosomal protein L37-A OS=Saccharomyces cerevisiae
"P51402" => "L37X", # 60S ribosomal protein L37-B OS=Saccharomyces cerevisiae
"P49167" => "L38", # 60S ribosomal protein L38 OS=Saccharomyces cerevisiae
"P04650" => "L39", # 60S ribosomal protein L39 OS=Saccharomyces cerevisiae
"P14126" => "L3", # 60S ribosomal protein L3 OS=Saccharomyces cerevisiae
"P14796" => "L40", # 60S ribosomal protein L40 OS=Saccharomyces cerevisiae
"P05746" => "L41", # 60S ribosomal protein L41 OS=Saccharomyces cerevisiae
"P49631" => "L43", # 60S ribosomal protein L43 OS=Saccharomyces cerevisiae
"P02405" => "L42", # 60S ribosomal protein L42 OS=Saccharomyces cerevisiae
"P10664" => "L4X", # 60S ribosomal protein L4-A OS=Saccharomyces cerevisiae
"P49626" => "L4X", # 60S ribosomal protein L4-B OS=Saccharomyces cerevisiae
"P26321" => "L5", # 60S ribosomal protein L5 OS=Saccharomyces cerevisiae
"Q02326" => "L6X", # 60S ribosomal protein L6-A OS=Saccharomyces cerevisiae
"P05739" => "L6X", # 60S ribosomal protein L6-B OS=Saccharomyces cerevisiae
"P05737" => "L7X", # 60S ribosomal protein L7-A OS=Saccharomyces cerevisiae
"Q12213" => "L7X", # 60S ribosomal protein L7-B OS=Saccharomyces cerevisiae
"P17076" => "L8X", # 60S ribosomal protein L8-A OS=Saccharomyces cerevisiae
"P29453" => "L8X", # 60S ribosomal protein L8-B OS=Saccharomyces cerevisiae
"P05738" => "L9X", # 60S ribosomal protein L9-A OS=Saccharomyces cerevisiae
"P51401" => "L9X", # 60S ribosomal protein L9-B OS=Saccharomyces cerevisiae
"P05317" => "P0", # 60S acidic ribosomal protein P0 OS=Saccharomyces cerevisiae
"P05318" => "P1X", # 60S acidic ribosomal protein P1-alpha OS=Saccharomyces cerevisiae
"P05319" => "P2X", # 60S acidic ribosomal protein P2-alpha OS=Saccharomyces cerevisiae
"P10622" => "P1X", # 60S acidic ribosomal protein P1-beta OS=Saccharomyces cerevisiae
"P02400" => "P2X", # 60S acidic ribosomal protein P2-beta OS=Saccharomyces cerevisiae


);

%acc2pro_human40S = (
# Human 40S
"P15880" => "RS2", # >sp|P15880|RS2_HUMAN 40S ribosomal protein S2 OS=Homo sapiens GN=RPS2 PE=1 SV=2
"P23396" => "RS3", # >sp|P23396|RS3_HUMAN 40S ribosomal protein S3 OS=Homo sapiens GN=RPS3 PE=1 SV=2
"P61247" => "RS3A", # >sp|P61247|RS3A_HUMAN 40S ribosomal protein S3a OS=Homo sapiens GN=RPS3A PE=1 SV=2
"P62701" => "RS4X", # >sp|P62701|RS4X_HUMAN 40S ribosomal protein S4, X isoform OS=Homo sapiens GN=RPS4X PE=1 SV=2
"P22090" => "RS4Y1", # >sp|P22090|RS4Y1_HUMAN 40S ribosomal protein S4, Y isoform 1 OS=Homo sapiens GN=RPS4Y1 PE=2 SV=2
"Q8TD47" => "RS4Y2", # >sp|Q8TD47|RS4Y2_HUMAN 40S ribosomal protein S4, Y isoform 2 OS=Homo sapiens GN=RPS4Y2 PE=1 SV=3
"P46782" => "RS5", # >sp|P46782|RS5_HUMAN 40S ribosomal protein S5 OS=Homo sapiens GN=RPS5 PE=1 SV=4
"P62753" => "RS6", # >sp|P62753|RS6_HUMAN 40S ribosomal protein S6 OS=Homo sapiens GN=RPS6 PE=1 SV=1
"P62081" => "RS7", # >sp|P62081|RS7_HUMAN 40S ribosomal protein S7 OS=Homo sapiens GN=RPS7 PE=1 SV=1
"P62241" => "RS8", # >sp|P62241|RS8_HUMAN 40S ribosomal protein S8 OS=Homo sapiens GN=RPS8 PE=1 SV=2
"P46781" => "RS9", # >sp|P46781|RS9_HUMAN 40S ribosomal protein S9 OS=Homo sapiens GN=RPS9 PE=1 SV=3
"P46783" => "RS10", # >sp|P46783|RS10_HUMAN 40S ribosomal protein S10 OS=Homo sapiens GN=RPS10 PE=1 SV=1
"Q9NQ39" => "RS10L", # >sp|Q9NQ39|RS10L_HUMAN Putative 40S ribosomal protein S10-like OS=Homo sapiens GN=RPS10P5 PE=5 SV=1
"P62280" => "RS11", # >sp|P62280|RS11_HUMAN 40S ribosomal protein S11 OS=Homo sapiens GN=RPS11 PE=1 SV=3
"P25398" => "RS12", # >sp|P25398|RS12_HUMAN 40S ribosomal protein S12 OS=Homo sapiens GN=RPS12 PE=1 SV=3
"P62277" => "RS13", # >sp|P62277|RS13_HUMAN 40S ribosomal protein S13 OS=Homo sapiens GN=RPS13 PE=1 SV=2
"P62263" => "RS14", # >sp|P62263|RS14_HUMAN 40S ribosomal protein S14 OS=Homo sapiens GN=RPS14 PE=1 SV=3
"P62841" => "RS15", # >sp|P62841|RS15_HUMAN 40S ribosomal protein S15 OS=Homo sapiens GN=RPS15 PE=1 SV=2
"P62244" => "RS15A", # >sp|P62244|RS15A_HUMAN 40S ribosomal protein S15a OS=Homo sapiens GN=RPS15A PE=1 SV=2
"P62249" => "RS16", # >sp|P62249|RS16_HUMAN 40S ribosomal protein S16 OS=Homo sapiens GN=RPS16 PE=1 SV=2
"P08708" => "RS17", # >sp|P08708|RS17_HUMAN 40S ribosomal protein S17 OS=Homo sapiens GN=RPS17 PE=1 SV=2
"P62269" => "RS18", # >sp|P62269|RS18_HUMAN 40S ribosomal protein S18 OS=Homo sapiens GN=RPS18 PE=1 SV=3
"P39019" => "RS19", # >sp|P39019|RS19_HUMAN 40S ribosomal protein S19 OS=Homo sapiens GN=RPS19 PE=1 SV=2
"P60866" => "RS20", # >sp|P60866|RS20_HUMAN 40S ribosomal protein S20 OS=Homo sapiens GN=RPS20 PE=1 SV=1
"P63220" => "RS21", # >sp|P63220|RS21_HUMAN 40S ribosomal protein S21 OS=Homo sapiens GN=RPS21 PE=1 SV=1
"P62266" => "RS23", # >sp|P62266|RS23_HUMAN 40S ribosomal protein S23 OS=Homo sapiens GN=RPS23 PE=1 SV=3
"P62847" => "RS24", # >sp|P62847|RS24_HUMAN 40S ribosomal protein S24 OS=Homo sapiens GN=RPS24 PE=1 SV=1
"P62851" => "RS25", # >sp|P62851|RS25_HUMAN 40S ribosomal protein S25 OS=Homo sapiens GN=RPS25 PE=1 SV=1
"P62854" => "RS26", # >sp|P62854|RS26_HUMAN 40S ribosomal protein S26 OS=Homo sapiens GN=RPS26 PE=1 SV=3
"Q5JNZ5" => "RS26L", # >sp|Q5JNZ5|RS26L_HUMAN Putative 40S ribosomal protein S26-like 1 OS=Homo sapiens GN=RPS26P11 PE=5 SV=1
"P42677" => "RS27", # >sp|P42677|RS27_HUMAN 40S ribosomal protein S27 OS=Homo sapiens GN=RPS27 PE=1 SV=3
"P62979" => "RS27A", # >sp|P62979|RS27A_HUMAN 40S ribosomal protein S27a OS=Homo sapiens GN=RPS27A PE=1 SV=1
"Q71UM5" => "RS27L", # >sp|Q71UM5|RS27L_HUMAN 40S ribosomal protein S27-like OS=Homo sapiens GN=RPS27L PE=1 SV=3
"P62857" => "RS28", # >sp|P62857|RS28_HUMAN 40S ribosomal protein S28 OS=Homo sapiens GN=RPS28 PE=1 SV=1
"P62273" => "RS29", # >sp|P62273|RS29_HUMAN 40S ribosomal protein S29 OS=Homo sapiens GN=RPS29 PE=1 SV=2
"P62861" => "RS30", # >sp|P62861|RS30_HUMAN 40S ribosomal protein S30 OS=Homo sapiens GN=FAU PE=1 SV=1
"P08865" => "RSSA", # >sp|P08865|RSSA_HUMAN 40S ribosomal protein SA OS=Homo sapiens GN=RPSA PE=1 SV=4
);

%acc2pro_human60S = (
# Human 60S
"P39023" => "RL3", # >sp|P39023|RL3_HUMAN 60S ribosomal protein L3 OS=Homo sapiens GN=RPL3 PE=1 SV=2
"Q92901" => "RL3L", # >sp|Q92901|RL3L_HUMAN 60S ribosomal protein L3-like OS=Homo sapiens GN=RPL3L PE=1 SV=3
"P36578" => "RL4", # >sp|P36578|RL4_HUMAN 60S ribosomal protein L4 OS=Homo sapiens GN=RPL4 PE=1 SV=5
"P46777" => "RL5", # >sp|P46777|RL5_HUMAN 60S ribosomal protein L5 OS=Homo sapiens GN=RPL5 PE=1 SV=3
"Q02878" => "RL6", # >sp|Q02878|RL6_HUMAN 60S ribosomal protein L6 OS=Homo sapiens GN=RPL6 PE=1 SV=3
"P18124" => "RL7", # >sp|P18124|RL7_HUMAN 60S ribosomal protein L7 OS=Homo sapiens GN=RPL7 PE=1 SV=1
"P62424" => "RL7A", # >sp|P62424|RL7A_HUMAN 60S ribosomal protein L7a OS=Homo sapiens GN=RPL7A PE=1 SV=2
"Q6DKI1" => "RL7L", # >sp|Q6DKI1|RL7L_HUMAN 60S ribosomal protein L7-like 1 OS=Homo sapiens GN=RPL7L1 PE=1 SV=1
"P62917" => "RL8", # >sp|P62917|RL8_HUMAN 60S ribosomal protein L8 OS=Homo sapiens GN=RPL8 PE=1 SV=2
"P32969" => "RL9", # >sp|P32969|RL9_HUMAN 60S ribosomal protein L9 OS=Homo sapiens GN=RPL9 PE=1 SV=1
"P27635" => "RL10", # >sp|P27635|RL10_HUMAN 60S ribosomal protein L10 OS=Homo sapiens GN=RPL10 PE=1 SV=4
"P62906" => "RL10A", # >sp|P62906|RL10A_HUMAN 60S ribosomal protein L10a OS=Homo sapiens GN=RPL10A PE=1 SV=2
"Q96L21" => "RL10L", # >sp|Q96L21|RL10L_HUMAN 60S ribosomal protein L10-like OS=Homo sapiens GN=RPL10L PE=1 SV=3
"P62913" => "RL11", # >sp|P62913|RL11_HUMAN 60S ribosomal protein L11 OS=Homo sapiens GN=RPL11 PE=1 SV=2
"P30050" => "RL12", # >sp|P30050|RL12_HUMAN 60S ribosomal protein L12 OS=Homo sapiens GN=RPL12 PE=1 SV=1
"P26373" => "RL13", # >sp|P26373|RL13_HUMAN 60S ribosomal protein L13 OS=Homo sapiens GN=RPL13 PE=1 SV=4
"P40429" => "RL13A", # >sp|P40429|RL13A_HUMAN 60S ribosomal protein L13a OS=Homo sapiens GN=RPL13A PE=1 SV=2
#"Q6NVV1" => "R13AX", # >sp|Q6NVV1|R13AX_HUMAN Putative 60S ribosomal protein L13a-like MGC87657 OS=Homo sapiens PE=5 SV=1
"Q6NVV1" => "RL13AX", # >sp|Q6NVV1|R13AX_HUMAN Putative 60S ribosomal protein L13a-like MGC87657 OS=Homo sapiens PE=5 SV=1
"P50914" => "RL14", # >sp|P50914|RL14_HUMAN 60S ribosomal protein L14 OS=Homo sapiens GN=RPL14 PE=1 SV=4
"P61313" => "RL15", # >sp|P61313|RL15_HUMAN 60S ribosomal protein L15 OS=Homo sapiens GN=RPL15 PE=1 SV=2
"P18621" => "RL17", # >sp|P18621|RL17_HUMAN 60S ribosomal protein L17 OS=Homo sapiens GN=RPL17 PE=1 SV=3
"Q07020" => "RL18", # >sp|Q07020|RL18_HUMAN 60S ribosomal protein L18 OS=Homo sapiens GN=RPL18 PE=1 SV=2
"Q02543" => "RL18A", # >sp|Q02543|RL18A_HUMAN 60S ribosomal protein L18a OS=Homo sapiens GN=RPL18A PE=1 SV=2
"P84098" => "RL19", # >sp|P84098|RL19_HUMAN 60S ribosomal protein L19 OS=Homo sapiens GN=RPL19 PE=1 SV=1
"P46778" => "RL21", # >sp|P46778|RL21_HUMAN 60S ribosomal protein L21 OS=Homo sapiens GN=RPL21 PE=1 SV=2
"P35268" => "RL22", # >sp|P35268|RL22_HUMAN 60S ribosomal protein L22 OS=Homo sapiens GN=RPL22 PE=1 SV=2
"Q6P5R6" => "RL22L", # >sp|Q6P5R6|RL22L_HUMAN 60S ribosomal protein L22-like 1 OS=Homo sapiens GN=RPL22L1 PE=1 SV=2
"P62829" => "RL23", # >sp|P62829|RL23_HUMAN 60S ribosomal protein L23 OS=Homo sapiens GN=RPL23 PE=1 SV=1
"P62750" => "RL23A", # >sp|P62750|RL23A_HUMAN 60S ribosomal protein L23a OS=Homo sapiens GN=RPL23A PE=1 SV=1
"P83731" => "RL24", # >sp|P83731|RL24_HUMAN 60S ribosomal protein L24 OS=Homo sapiens GN=RPL24 PE=1 SV=1
"P61254" => "RL26", # >sp|P61254|RL26_HUMAN 60S ribosomal protein L26 OS=Homo sapiens GN=RPL26 PE=1 SV=1
"Q9UNX3" => "RL26L", # >sp|Q9UNX3|RL26L_HUMAN 60S ribosomal protein L26-like 1 OS=Homo sapiens GN=RPL26L1 PE=1 SV=1
"P61353" => "RL27", # >sp|P61353|RL27_HUMAN 60S ribosomal protein L27 OS=Homo sapiens GN=RPL27 PE=1 SV=2
"P46776" => "RL27A", # >sp|P46776|RL27A_HUMAN 60S ribosomal protein L27a OS=Homo sapiens GN=RPL27A PE=1 SV=2
"P46779" => "RL28", # >sp|P46779|RL28_HUMAN 60S ribosomal protein L28 OS=Homo sapiens GN=RPL28 PE=1 SV=3
"P47914" => "RL29", # >sp|P47914|RL29_HUMAN 60S ribosomal protein L29 OS=Homo sapiens GN=RPL29 PE=1 SV=2
"P62888" => "RL30", # >sp|P62888|RL30_HUMAN 60S ribosomal protein L30 OS=Homo sapiens GN=RPL30 PE=1 SV=2
"P62899" => "RL31", # >sp|P62899|RL31_HUMAN 60S ribosomal protein L31 OS=Homo sapiens GN=RPL31 PE=1 SV=1
"P62910" => "RL32", # >sp|P62910|RL32_HUMAN 60S ribosomal protein L32 OS=Homo sapiens GN=RPL32 PE=1 SV=2
"P49207" => "RL34", # >sp|P49207|RL34_HUMAN 60S ribosomal protein L34 OS=Homo sapiens GN=RPL34 PE=1 SV=3
"P42766" => "RL35", # >sp|P42766|RL35_HUMAN 60S ribosomal protein L35 OS=Homo sapiens GN=RPL35 PE=1 SV=2
"P18077" => "RL35A", # >sp|P18077|RL35A_HUMAN 60S ribosomal protein L35a OS=Homo sapiens GN=RPL35A PE=1 SV=2
"Q9Y3U8" => "RL36", # >sp|Q9Y3U8|RL36_HUMAN 60S ribosomal protein L36 OS=Homo sapiens GN=RPL36 PE=1 SV=3
"P83881" => "RL36A", # >sp|P83881|RL36A_HUMAN 60S ribosomal protein L36a OS=Homo sapiens GN=RPL36A PE=1 SV=2
"Q969Q0" => "RL36L", # >sp|Q969Q0|RL36L_HUMAN 60S ribosomal protein L36a-like OS=Homo sapiens GN=RPL36AL PE=1 SV=3
"P0C6E6" => "RL36X", # >sp|P0C6E6|RL36X_HUMAN Putative 60S ribosomal protein L36-like 1 OS=Homo sapiens PE=5 SV=1
"P61927" => "RL37", # >sp|P61927|RL37_HUMAN 60S ribosomal protein L37 OS=Homo sapiens GN=RPL37 PE=1 SV=2
"P61513" => "RL37A", # >sp|P61513|RL37A_HUMAN 60S ribosomal protein L37a OS=Homo sapiens GN=RPL37A PE=1 SV=2
"A6NKH3" => "RL37L", # >sp|A6NKH3|RL37L_HUMAN Putative 60S ribosomal protein L37a-like OS=Homo sapiens GN=RPL37L PE=5 SV=2
"P63173" => "RL38", # >sp|P63173|RL38_HUMAN 60S ribosomal protein L38 OS=Homo sapiens GN=RPL38 PE=1 SV=2
"P62891" => "RL39", # >sp|P62891|RL39_HUMAN 60S ribosomal protein L39 OS=Homo sapiens GN=RPL39 PE=2 SV=2
"Q96EH5" => "RL39L", # >sp|Q96EH5|RL39L_HUMAN 60S ribosomal protein L39-like OS=Homo sapiens GN=RPL39L PE=1 SV=3
#"Q59GN2" => "R39L5", # >sp|Q59GN2|R39L5_HUMAN Putative 60S ribosomal protein L39-like 5 OS=Homo sapiens GN=RPL39P5 PE=5 SV=2
"Q59GN2" => "RL39L5", # >sp|Q59GN2|R39L5_HUMAN Putative 60S ribosomal protein L39-like 5 OS=Homo sapiens GN=RPL39P5 PE=5 SV=2
"P62987" => "RL40", # >sp|P62987|RL40_HUMAN 60S ribosomal protein L40 OS=Homo sapiens GN=UBA52 PE=1 SV=1
"P62945" => "RL41", # >sp|P62945|RL41_HUMAN 60S ribosomal protein L41 OS=Homo sapiens GN=RPL41 PE=2 SV=1
"P05388" => "RLA0", # >sp|P05388|RLA0_HUMAN 60S acidic ribosomal protein P0 OS=Homo sapiens GN=RPLP0 PE=1 SV=1
"Q8NHW5" => "RLA0L", # >sp|Q8NHW5|RLA0L_HUMAN 60S acidic ribosomal protein P0-like OS=Homo sapiens PE=1 SV=1
"P05386" => "RLA1", # >sp|P05386|RLA1_HUMAN 60S acidic ribosomal protein P1 OS=Homo sapiens GN=RPLP1 PE=1 SV=1
"P05387" => "RLA2", # >sp|P05387|RLA2_HUMAN 60S acidic ribosomal protein P2 OS=Homo sapiens GN=RPLP2 PE=1 SV=1
);

%acc2pro = (
	%acc2pro_ecoli30S,
	%acc2pro_ecoli50S,
	%acc2pro_yeast40S,
	%acc2pro_yeast60S,
	%acc2pro_human40S,
	%acc2pro_human60S
);


%acc2pro_all = (
	%acc2pro,

# Human Ribosomal Kinase and Associate Proteins
"Q15418" => "KS6A1", # >sp|Q15418|KS6A1_HUMAN Ribosomal protein S6 kinase alpha-1 OS=Homo sapiens GN=RPS6KA1 PE=1 SV=2
"Q15349" => "KS6A2", # >sp|Q15349|KS6A2_HUMAN Ribosomal protein S6 kinase alpha-2 OS=Homo sapiens GN=RPS6KA2 PE=1 SV=2
"P51812" => "KS6A3", # >sp|P51812|KS6A3_HUMAN Ribosomal protein S6 kinase alpha-3 OS=Homo sapiens GN=RPS6KA3 PE=1 SV=1
"O75676" => "KS6A4", # >sp|O75676|KS6A4_HUMAN Ribosomal protein S6 kinase alpha-4 OS=Homo sapiens GN=RPS6KA4 PE=1 SV=1
"O75582" => "KS6A5", # >sp|O75582|KS6A5_HUMAN Ribosomal protein S6 kinase alpha-5 OS=Homo sapiens GN=RPS6KA5 PE=1 SV=1
"Q9UK32" => "KS6A6", # >sp|Q9UK32|KS6A6_HUMAN Ribosomal protein S6 kinase alpha-6 OS=Homo sapiens GN=RPS6KA6 PE=1 SV=1
"P23443" => "KS6B1", # >sp|P23443|KS6B1_HUMAN Ribosomal protein S6 kinase beta-1 OS=Homo sapiens GN=RPS6KB1 PE=1 SV=2
"Q9UBS0" => "KS6B2", # >sp|Q9UBS0|KS6B2_HUMAN Ribosomal protein S6 kinase beta-2 OS=Homo sapiens GN=RPS6KB2 PE=1 SV=2
"Q96S38" => "KS6C1", # >sp|Q96S38|KS6C1_HUMAN Ribosomal protein S6 kinase delta-1 OS=Homo sapiens GN=RPS6KC1 PE=1 SV=2
"Q96D46" => "NMD3", # >sp|Q96D46|NMD3_HUMAN 60S ribosomal export protein NMD3 OS=Homo sapiens GN=NMD3 PE=1 SV=1
"P46087" => "NOP2", # >sp|P46087|NOP2_HUMAN Putative ribosomal RNA methyltransferase NOP2 OS=Homo sapiens GN=NOP2 PE=1 SV=2
"O60287" => "NPA1P", # >sp|O60287|NPA1P_HUMAN Nucleolar pre-ribosomal-associated protein 1 OS=Homo sapiens GN=URB1 PE=1 SV=4
"Q8IXN7" => "RIMKA", # >sp|Q8IXN7|RIMKA_HUMAN Ribosomal protein S6 modification-like protein A OS=Homo sapiens GN=RIMKLA PE=2 SV=2
"Q9ULI2" => "RIMKB", # >sp|Q9ULI2|RIMKB_HUMAN Ribosomal protein S6 modification-like protein B OS=Homo sapiens GN=RIMKLB PE=2 SV=2
"Q9Y6S9" => "RPKL1", # >sp|Q9Y6S9|RPKL1_HUMAN Ribosomal protein S6 kinase-like 1 OS=Homo sapiens GN=RPS6KL1 PE=2 SV=1
"Q9UET6" => "RRMJ1", # >sp|Q9UET6|RRMJ1_HUMAN Putative ribosomal RNA methyltransferase 1 OS=Homo sapiens GN=FTSJ1 PE=1 SV=2
"Q9UI43" => "RRMJ2", # >sp|Q9UI43|RRMJ2_HUMAN Putative ribosomal RNA methyltransferase 2 OS=Homo sapiens GN=FTSJ2 PE=1 SV=1
"Q14684" => "RRP1B", # >sp|Q14684|RRP1B_HUMAN Ribosomal RNA processing protein 1 homolog B OS=Homo sapiens GN=RRP1B PE=1 SV=3
"P56182" => "RRP1", # >sp|P56182|RRP1_HUMAN Ribosomal RNA processing protein 1 homolog A OS=Homo sapiens GN=RRP1 PE=1 SV=1
"Q9Y3A4" => "RRP7A", # >sp|Q9Y3A4|RRP7A_HUMAN Ribosomal RNA-processing protein 7 homolog A OS=Homo sapiens GN=RRP7A PE=1 SV=2
"Q9NSQ0" => "RRP7B", # >sp|Q9NSQ0|RRP7B_HUMAN Putative ribosomal RNA-processing protein 7 homolog B OS=Homo sapiens GN=RRP7B PE=5 SV=1
"O43159" => "RRP8", # >sp|O43159|RRP8_HUMAN Ribosomal RNA-processing protein 8 OS=Homo sapiens GN=RRP8 PE=1 SV=2

# Human Mitochondrial Ribosome 28S
"Q9Y399" => "RT02", # >sp|Q9Y399|RT02_HUMAN 28S ribosomal protein S2, mitochondrial OS=Homo sapiens GN=MRPS2 PE=1 SV=1
"P82675" => "RT05", # >sp|P82675|RT05_HUMAN 28S ribosomal protein S5, mitochondrial OS=Homo sapiens GN=MRPS5 PE=1 SV=2
"P82932" => "RT06", # >sp|P82932|RT06_HUMAN 28S ribosomal protein S6, mitochondrial OS=Homo sapiens GN=MRPS6 PE=1 SV=3
"Q9Y2R9" => "RT07", # >sp|Q9Y2R9|RT07_HUMAN 28S ribosomal protein S7, mitochondrial OS=Homo sapiens GN=MRPS7 PE=1 SV=2
"P82933" => "RT09", # >sp|P82933|RT09_HUMAN 28S ribosomal protein S9, mitochondrial OS=Homo sapiens GN=MRPS9 PE=1 SV=2
"P82664" => "RT10", # >sp|P82664|RT10_HUMAN 28S ribosomal protein S10, mitochondrial OS=Homo sapiens GN=MRPS10 PE=1 SV=2
"P82912" => "RT11", # >sp|P82912|RT11_HUMAN 28S ribosomal protein S11, mitochondrial OS=Homo sapiens GN=MRPS11 PE=1 SV=2
"O15235" => "RT12", # >sp|O15235|RT12_HUMAN 28S ribosomal protein S12, mitochondrial OS=Homo sapiens GN=MRPS12 PE=1 SV=1
"O60783" => "RT14", # >sp|O60783|RT14_HUMAN 28S ribosomal protein S14, mitochondrial OS=Homo sapiens GN=MRPS14 PE=1 SV=1
"P82914" => "RT15", # >sp|P82914|RT15_HUMAN 28S ribosomal protein S15, mitochondrial OS=Homo sapiens GN=MRPS15 PE=1 SV=1
"Q9Y3D3" => "RT16", # >sp|Q9Y3D3|RT16_HUMAN 28S ribosomal protein S16, mitochondrial OS=Homo sapiens GN=MRPS16 PE=1 SV=1
"Q9Y2R5" => "RT17", # >sp|Q9Y2R5|RT17_HUMAN 28S ribosomal protein S17, mitochondrial OS=Homo sapiens GN=MRPS17 PE=1 SV=1
"Q9NVS2" => "RT18A", # >sp|Q9NVS2|RT18A_HUMAN 28S ribosomal protein S18a, mitochondrial OS=Homo sapiens GN=MRPS18A PE=1 SV=1
"Q9Y676" => "RT18B", # >sp|Q9Y676|RT18B_HUMAN 28S ribosomal protein S18b, mitochondrial OS=Homo sapiens GN=MRPS18B PE=1 SV=1
"Q9Y3D5" => "RT18C", # >sp|Q9Y3D5|RT18C_HUMAN 28S ribosomal protein S18c, mitochondrial OS=Homo sapiens GN=MRPS18C PE=1 SV=1
"P82921" => "RT21", # >sp|P82921|RT21_HUMAN 28S ribosomal protein S21, mitochondrial OS=Homo sapiens GN=MRPS21 PE=1 SV=2
"P82650" => "RT22", # >sp|P82650|RT22_HUMAN 28S ribosomal protein S22, mitochondrial OS=Homo sapiens GN=MRPS22 PE=1 SV=1
"Q9Y3D9" => "RT23", # >sp|Q9Y3D9|RT23_HUMAN 28S ribosomal protein S23, mitochondrial OS=Homo sapiens GN=MRPS23 PE=1 SV=2
"Q96EL2" => "RT24", # >sp|Q96EL2|RT24_HUMAN 28S ribosomal protein S24, mitochondrial OS=Homo sapiens GN=MRPS24 PE=1 SV=1
"P82663" => "RT25", # >sp|P82663|RT25_HUMAN 28S ribosomal protein S25, mitochondrial OS=Homo sapiens GN=MRPS25 PE=1 SV=1
"Q9BYN8" => "RT26", # >sp|Q9BYN8|RT26_HUMAN 28S ribosomal protein S26, mitochondrial OS=Homo sapiens GN=MRPS26 PE=1 SV=1
"Q92552" => "RT27", # >sp|Q92552|RT27_HUMAN 28S ribosomal protein S27, mitochondrial OS=Homo sapiens GN=MRPS27 PE=1 SV=3
"Q9Y2Q9" => "RT28", # >sp|Q9Y2Q9|RT28_HUMAN 28S ribosomal protein S28, mitochondrial OS=Homo sapiens GN=MRPS28 PE=1 SV=1
"P51398" => "RT29", # >sp|P51398|RT29_HUMAN 28S ribosomal protein S29, mitochondrial OS=Homo sapiens GN=DAP3 PE=1 SV=1
"Q9NP92" => "RT30", # >sp|Q9NP92|RT30_HUMAN 28S ribosomal protein S30, mitochondrial OS=Homo sapiens GN=MRPS30 PE=1 SV=2
"Q92665" => "RT31", # >sp|Q92665|RT31_HUMAN 28S ribosomal protein S31, mitochondrial OS=Homo sapiens GN=MRPS31 PE=1 SV=3
"Q9Y291" => "RT33", # >sp|Q9Y291|RT33_HUMAN 28S ribosomal protein S33, mitochondrial OS=Homo sapiens GN=MRPS33 PE=1 SV=1
"P82930" => "RT34", # >sp|P82930|RT34_HUMAN 28S ribosomal protein S34, mitochondrial OS=Homo sapiens GN=MRPS34 PE=1 SV=2
"P82673" => "RT35", # >sp|P82673|RT35_HUMAN 28S ribosomal protein S35, mitochondrial OS=Homo sapiens GN=MRPS35 PE=1 SV=1
"P82909" => "RT36", # >sp|P82909|RT36_HUMAN 28S ribosomal protein S36, mitochondrial OS=Homo sapiens GN=MRPS36 PE=1 SV=2

# Human Mitochondrial Ribosome 39S
"Q9BYD6" => "RM01", # >sp|Q9BYD6|RM01_HUMAN 39S ribosomal protein L1, mitochondrial OS=Homo sapiens GN=MRPL1 PE=1 SV=2
"Q5T653" => "RM02", # >sp|Q5T653|RM02_HUMAN 39S ribosomal protein L2, mitochondrial OS=Homo sapiens GN=MRPL2 PE=1 SV=2
"P09001" => "RM03", # >sp|P09001|RM03_HUMAN 39S ribosomal protein L3, mitochondrial OS=Homo sapiens GN=MRPL3 PE=1 SV=1
"Q9BYD3" => "RM04", # >sp|Q9BYD3|RM04_HUMAN 39S ribosomal protein L4, mitochondrial OS=Homo sapiens GN=MRPL4 PE=1 SV=1
"Q9BYD2" => "RM09", # >sp|Q9BYD2|RM09_HUMAN 39S ribosomal protein L9, mitochondrial OS=Homo sapiens GN=MRPL9 PE=1 SV=2
"Q7Z7H8" => "RM10", # >sp|Q7Z7H8|RM10_HUMAN 39S ribosomal protein L10, mitochondrial OS=Homo sapiens GN=MRPL10 PE=1 SV=3
"Q9Y3B7" => "RM11", # >sp|Q9Y3B7|RM11_HUMAN 39S ribosomal protein L11, mitochondrial OS=Homo sapiens GN=MRPL11 PE=1 SV=1
"P52815" => "RM12", # >sp|P52815|RM12_HUMAN 39S ribosomal protein L12, mitochondrial OS=Homo sapiens GN=MRPL12 PE=1 SV=2
"Q9BYD1" => "RM13", # >sp|Q9BYD1|RM13_HUMAN 39S ribosomal protein L13, mitochondrial OS=Homo sapiens GN=MRPL13 PE=1 SV=1
"Q6P1L8" => "RM14", # >sp|Q6P1L8|RM14_HUMAN 39S ribosomal protein L14, mitochondrial OS=Homo sapiens GN=MRPL14 PE=1 SV=1
"Q9P015" => "RM15", # >sp|Q9P015|RM15_HUMAN 39S ribosomal protein L15, mitochondrial OS=Homo sapiens GN=MRPL15 PE=1 SV=1
"Q9NX20" => "RM16", # >sp|Q9NX20|RM16_HUMAN 39S ribosomal protein L16, mitochondrial OS=Homo sapiens GN=MRPL16 PE=1 SV=1
"Q9NRX2" => "RM17", # >sp|Q9NRX2|RM17_HUMAN 39S ribosomal protein L17, mitochondrial OS=Homo sapiens GN=MRPL17 PE=1 SV=1
"Q9H0U6" => "RM18", # >sp|Q9H0U6|RM18_HUMAN 39S ribosomal protein L18, mitochondrial OS=Homo sapiens GN=MRPL18 PE=1 SV=1
"P49406" => "RM19", # >sp|P49406|RM19_HUMAN 39S ribosomal protein L19, mitochondrial OS=Homo sapiens GN=MRPL19 PE=1 SV=2
"Q9BYC9" => "RM20", # >sp|Q9BYC9|RM20_HUMAN 39S ribosomal protein L20, mitochondrial OS=Homo sapiens GN=MRPL20 PE=1 SV=1
"Q7Z2W9" => "RM21", # >sp|Q7Z2W9|RM21_HUMAN 39S ribosomal protein L21, mitochondrial OS=Homo sapiens GN=MRPL21 PE=1 SV=2
"Q9NWU5" => "RM22", # >sp|Q9NWU5|RM22_HUMAN 39S ribosomal protein L22, mitochondrial OS=Homo sapiens GN=MRPL22 PE=1 SV=1
"Q16540" => "RM23", # >sp|Q16540|RM23_HUMAN 39S ribosomal protein L23, mitochondrial OS=Homo sapiens GN=MRPL23 PE=1 SV=1
"Q96A35" => "RM24", # >sp|Q96A35|RM24_HUMAN 39S ribosomal protein L24, mitochondrial OS=Homo sapiens GN=MRPL24 PE=1 SV=1
"Q9P0M9" => "RM27", # >sp|Q9P0M9|RM27_HUMAN 39S ribosomal protein L27, mitochondrial OS=Homo sapiens GN=MRPL27 PE=1 SV=1
"Q13084" => "RM28", # >sp|Q13084|RM28_HUMAN 39S ribosomal protein L28, mitochondrial OS=Homo sapiens GN=MRPL28 PE=1 SV=4
"Q8TCC3" => "RM30", # >sp|Q8TCC3|RM30_HUMAN 39S ribosomal protein L30, mitochondrial OS=Homo sapiens GN=MRPL30 PE=1 SV=1
"Q9BYC8" => "RM32", # >sp|Q9BYC8|RM32_HUMAN 39S ribosomal protein L32, mitochondrial OS=Homo sapiens GN=MRPL32 PE=1 SV=1
"O75394" => "RM33", # >sp|O75394|RM33_HUMAN 39S ribosomal protein L33, mitochondrial OS=Homo sapiens GN=MRPL33 PE=1 SV=1
"Q9BQ48" => "RM34", # >sp|Q9BQ48|RM34_HUMAN 39S ribosomal protein L34, mitochondrial OS=Homo sapiens GN=MRPL34 PE=1 SV=1
"Q9NZE8" => "RM35", # >sp|Q9NZE8|RM35_HUMAN 39S ribosomal protein L35, mitochondrial OS=Homo sapiens GN=MRPL35 PE=2 SV=3
"Q9P0J6" => "RM36", # >sp|Q9P0J6|RM36_HUMAN 39S ribosomal protein L36, mitochondrial OS=Homo sapiens GN=MRPL36 PE=2 SV=1
"Q9BZE1" => "RM37", # >sp|Q9BZE1|RM37_HUMAN 39S ribosomal protein L37, mitochondrial OS=Homo sapiens GN=MRPL37 PE=1 SV=2
"Q96DV4" => "RM38", # >sp|Q96DV4|RM38_HUMAN 39S ribosomal protein L38, mitochondrial OS=Homo sapiens GN=MRPL38 PE=1 SV=2
"Q9NYK5" => "RM39", # >sp|Q9NYK5|RM39_HUMAN 39S ribosomal protein L39, mitochondrial OS=Homo sapiens GN=MRPL39 PE=1 SV=3
"Q9NQ50" => "RM40", # >sp|Q9NQ50|RM40_HUMAN 39S ribosomal protein L40, mitochondrial OS=Homo sapiens GN=MRPL40 PE=1 SV=1
"Q8IXM3" => "RM41", # >sp|Q8IXM3|RM41_HUMAN 39S ribosomal protein L41, mitochondrial OS=Homo sapiens GN=MRPL41 PE=1 SV=1
"Q9Y6G3" => "RM42", # >sp|Q9Y6G3|RM42_HUMAN 39S ribosomal protein L42, mitochondrial OS=Homo sapiens GN=MRPL42 PE=1 SV=1
"Q8N983" => "RM43", # >sp|Q8N983|RM43_HUMAN 39S ribosomal protein L43, mitochondrial OS=Homo sapiens GN=MRPL43 PE=1 SV=1
"Q9H9J2" => "RM44", # >sp|Q9H9J2|RM44_HUMAN 39S ribosomal protein L44, mitochondrial OS=Homo sapiens GN=MRPL44 PE=1 SV=1
"Q9BRJ2" => "RM45", # >sp|Q9BRJ2|RM45_HUMAN 39S ribosomal protein L45, mitochondrial OS=Homo sapiens GN=MRPL45 PE=1 SV=2
"Q9H2W6" => "RM46", # >sp|Q9H2W6|RM46_HUMAN 39S ribosomal protein L46, mitochondrial OS=Homo sapiens GN=MRPL46 PE=1 SV=1
"Q9HD33" => "RM47", # >sp|Q9HD33|RM47_HUMAN 39S ribosomal protein L47, mitochondrial OS=Homo sapiens GN=MRPL47 PE=1 SV=2
"Q96GC5" => "RM48", # >sp|Q96GC5|RM48_HUMAN 39S ribosomal protein L48, mitochondrial OS=Homo sapiens GN=MRPL48 PE=1 SV=2
"Q13405" => "RM49", # >sp|Q13405|RM49_HUMAN 39S ribosomal protein L49, mitochondrial OS=Homo sapiens GN=MRPL49 PE=1 SV=1
"Q8N5N7" => "RM50", # >sp|Q8N5N7|RM50_HUMAN 39S ribosomal protein L50, mitochondrial OS=Homo sapiens GN=MRPL50 PE=1 SV=2
"Q4U2R6" => "RM51", # >sp|Q4U2R6|RM51_HUMAN 39S ribosomal protein L51, mitochondrial OS=Homo sapiens GN=MRPL51 PE=1 SV=1
"Q86TS9" => "RM52", # >sp|Q86TS9|RM52_HUMAN 39S ribosomal protein L52, mitochondrial OS=Homo sapiens GN=MRPL52 PE=2 SV=2
"Q96EL3" => "RM53", # >sp|Q96EL3|RM53_HUMAN 39S ribosomal protein L53, mitochondrial OS=Homo sapiens GN=MRPL53 PE=1 SV=1
"Q6P161" => "RM54", # >sp|Q6P161|RM54_HUMAN 39S ribosomal protein L54, mitochondrial OS=Homo sapiens GN=MRPL54 PE=1 SV=1
"Q7Z7F7" => "RM55", # >sp|Q7Z7F7|RM55_HUMAN 39S ribosomal protein L55, mitochondrial OS=Homo sapiens GN=MRPL55 PE=1 SV=1

);

# protein -> location

%pro2loc30Splot = (
	"XXX" => 0,
	"S2" => 1,
	"S3" => 2,
	"S4" => 3,
	"S5" => 4,
	"S6" => 5,
	"S7" => 6,
	"S8" => 7,
	"S9" => 8,
	"S10" => 9,
	"S11" => 10,
	"S12" => 11,
	"S13" => 12,
	"S14" => 13,
	"S15" => 14,
	"S16" => 15,
	"S17" => 16,
	"S18" => 17,
	"S19" => 18,
	"S20L26" => 19,
	"S20" => 19,
	"S21" => 20,
);

%pro2loc30S = (
	"XXX" => 0,
	"S1" => 1,
	"S2" => 2,
	"S3" => 3,
	"S4" => 4,
	"S5" => 5,
	"S6" => 6,
	"S7" => 7,
	"S8" => 8,
	"S9" => 9,
	"S10" => 10,
	"S11" => 11,
	"S12" => 12,
	"S13" => 13,
	"S14" => 14,
	"S15" => 15,
	"S16" => 16,
	"S17" => 17,
	"S18" => 18,
	"S19" => 19,
	"S20L26" => 20,
	"S20" => 20,
	"S21" => 21,
	"S22" => 22,
);

%pro2loc50S = (
	"XXX" => 0,
	"L1" => 1,
	"L2" => 2,
	"L3" => 3,
	"L4" => 4,
	"L5" => 5,
	"L6" => 6,
	"L7" => 7,
	"L7L12" => 7,
	"L8" => 8,
	"L9" => 9,
	"L10" => 10,
	"L11" => 11,
	"L12" => 12,
	"L13" => 13,
	"L14" => 14,
	"L15" => 15,
	"L16" => 16,
	"L17" => 17,
	"L18" => 18,
	"L19" => 19,
	"L20" => 20,
	"L21" => 21,
	"L22" => 22,
	"L23" => 23,
	"L24" => 24,
	"L25" => 25,
	"L26" => 26,
	"L27" => 27,
	"L28" => 28,
	"L29" => 29,
	"L30" => 30,
	"L31" => 31,
	"L32" => 32,
	"L33" => 33,
	"L34" => 34,
	"L35" => 35,
	"L36" => 36,
);

%pro2loc50Splot = (
	"XXX" => 0,
	"L1" => 1,
	"L2" => 2,
	"L3" => 3,
	"L4" => 4,
	"L5" => 5,
	"L6" => 6,
	"L7" => 7,
	"L7L12" => 7,
	"L9" => 8,
	"L10" => 9,
	"L11" => 10,
	"L13" => 11,
	"L14" => 12,
	"L15" => 13,
	"L16" => 14,
	"L17" => 15,
	"L18" => 16,
	"L19" => 17,
	"L20" => 18,
	"L21" => 19,
	"L22" => 20,
	"L23" => 21,
	"L24" => 22,
	"L25" => 23,
	"L26" => 24,
	"L27" => 25,
	"L28" => 26,
	"L29" => 27,
	"L30" => 28,
	"L31" => 29,
	"L32" => 30,
	"L33" => 31,
	"L34" => 32,
	"L35" => 33,
	"L36" => 34,
);

%pro2loc70S = (
	"XXX" => 0,
	"S1" => 1,
	"S2" => 2,
	"S3" => 3,
	"S4" => 4,
	"S5" => 5,
	"S6" => 6,
	"S7" => 7,
	"S8" => 8,
	"S9" => 9,
	"S10" => 10,
	"S11" => 11,
	"S12" => 12,
	"S13" => 13,
	"S14" => 14,
	"S15" => 15,
	"S16" => 16,
	"S17" => 17,
	"S18" => 18,
	"S19" => 19,
	"S20" => 20, "S20L26" => 20,
	"S21" => 21,
	"S22" => 22,
	
	"L1" => 26,
	"L2" => 27,
	"L3" => 28,
	"L4" => 29,
	"L5" => 30,
	"L6" => 31,
	"L7" => 32,
	"L7L12" => 32,
	"L8" => 33,
	"L9" => 34,
	"L10" => 35,
	"L11" => 36,
	"L12" => 37,
	"L13" => 38,
	"L14" => 39,
	"L15" => 40,
	"L16" => 41,
	"L17" => 42,
	"L18" => 43,
	"L19" => 44,
	"L20" => 45,
	"L21" => 46,
	"L22" => 47,
	"L23" => 48,
	"L24" => 49,
	"L25" => 50,
	"L26" => 51,
	"L27" => 52,
	"L28" => 53,
	"L29" => 54,
	"L30" => 55,
	"L31" => 56,
	"L32" => 57,
	"L33" => 58,
	"L34" => 59,
	"L35" => 60,
	"L36" => 61,
);

%pro2loc30Schlamy = (
	"XXX" => 0,
	"S1" => 1,
	"S2" => 2,
	"S2b" => 3,
	"S3" => 4,
	"S4" => 5,
	"S5" => 6,
	"S6" => 7,
	"S7" => 8,
	"S8" => 9,
	"S9" => 10,
	"S10" => 11,
	"S11" => 12,
	"S12" => 13,
	"S13" => 14,
	"S14" => 15,
	"S15" => 16,
	"S16" => 17,
	"S17" => 18,
	"S18" => 19,
	"S19" => 20,
	"S20" => 21,
	"S21" => 22,
	"PSRP1" => 23,
	"PSRP3" => 24,
	"PSRP4" => 25,
	"PSRP7full" => 26,
);

%pro2loc40Shuman = (
	"XXX" => 0,
	"SA" => 1,
	"S2" => 2,
	"S3" => 3,
	"S3A" => 4,
	"S4A" => 5,
	"S4B" => 5,
	"S4X" => 5,
	"S5" => 6,
	"S6" => 7,
	"S7" => 8,
	"S8" => 9,
	"S9" => 10,
	"S10" => 11,
	"S11" => 12,
	"S12" => 13,
	"S13" => 14,
	"S14" => 15,
	"S15" => 16,
	"S15A" => 17,
	"S16" => 18,
	"S17" => 19,
	"S18" => 20,
	"S19" => 21,
	"S20" => 22,
	"S21" => 23,
	"S23" => 24,
	"S24" => 25,
	"S25" => 26,
	"S26" => 27,
	"S27" => 28,
	"S27A" => 29,
	"S28" => 30,
	"S29" => 31,
	"S30" => 32,
);

%pro2loc60Shuman = (
	"XXX" => 0,
	"L3" => 1,
	"L4" => 2,
	"L5" => 3,
	"L6" => 4,
	"L7" => 5,
	"L7A" => 6,
	"L8" => 7,
	"L9" => 8,
	"L10" => 9,
	"L10A" => 10,
	"L11" => 11,
	"L12" => 12,
	"L13" => 13,
	"L13A" => 14,
	"L14" => 15,
	"L15" => 16,
	"L17" => 17,
	"L18" => 18,
	"L18A" => 19,
	"L19" => 20,
	"L21" => 21,
	"L22" => 22,
	"L23" => 23,
	"L23A" => 24,
	"L24" => 25,
	"L26" => 26,
	"L27" => 27,
	"L27A" => 28,
	"L28" => 29,
	"L29" => 30,
	"L30" => 31,
	"L31" => 32,
	"L32" => 33,
	"L34" => 34,
	"L35" => 35,
	"L35A" => 36,
	"L36" => 37,
	"L36A" => 38,
	"L37" => 39,
	"L37A" => 40,
	"L38" => 41,
	"L39" => 42,
	"L40" => 43,
	"L41" => 44,
	"LP0" => 45,
	"LP1" => 46,
	"LP2" => 47,
);

%pro2loc80Shuman = (
	"XXX" => 0,
	"SA" => 1,
	"S2" => 2,
	"S3" => 3,
	"S3A" => 4,
	"S4A" => 5,
	"S4B" => 5,
	"S4X" => 5,
	"S5" => 6,
	"S6" => 7,
	"S7" => 8,
	"S8" => 9,
	"S9" => 10,
	"S10" => 11,
	"S11" => 12,
	"S12" => 13,
	"S13" => 14,
	"S14" => 15,
	"S15" => 16,
	"S15A" => 17,
	"S16" => 18,
	"S17" => 19,
	"S18" => 20,
	"S19" => 21,
	"S20" => 22,
	"S21" => 23,
	"S23" => 24,
	"S24" => 25,
	"S25" => 26,
	"S26" => 27,
	"S27" => 28,
	"S27A" => 29,
	"S28" => 30,
	"S29" => 31,
	"S30" => 32,

	"L3" => 33,
	"L4" => 34,
	"L5" => 35,
	"L6" => 36,
	"L7" => 37,
	"L7A" => 38,
	"L8" => 39,
	"L9" => 40,
	"L10" => 41,
	"L10A" => 42,
	"L11" => 43,
	"L12" => 44,
	"L13" => 45,
	"L13A" => 46,
	"L14" => 47,
	"L15" => 48,
	"L17" => 49,
	"L18" => 50,
	"L18A" => 51,
	"L19" => 52,
	"L21" => 53,
	"L22" => 54,
	"L23" => 55,
	"L23A" => 56,
	"L24" => 57,
	"L26" => 58,
	"L27" => 59,
	"L27A" => 60,
	"L28" => 61,
	"L29" => 62,
	"L30" => 63,
	"L31" => 64,
	"L32" => 65,
	"L34" => 66,
	"L35" => 67,
	"L35A" => 68,
	"L36" => 69,
	"L36A" => 70,
	"L37" => 71,
	"L37A" => 72,
	"L38" => 73,
	"L39" => 74,
	"L40" => 75,
	"L41" => 76,
	"LP0" => 77,
	"LP1" => 78,
	"LP2" => 79,
);

%pro2loccof = (
        "XXX" => 0,
        "CsdA" => 1,
        "DbpA" => 2,
        "RrmJ" => 3,
        "SrmB" => 4,
);

%pro2loc40Syeast_old = (
	"XXX" => 0,
	"RS0A" => 1,
	"RS0B" => 2,
	"RS0X" => 3,
	"RS1A" => 4,
	"RS1B" => 5,
	"RS1X" => 6,
	"RS2" => 7,
	"RS3" => 8,
	"RS4A" => 9,
	"RS4B" => 10,
	"RS4X" => 11,
	"RS5" => 12,
	"RS6A" => 13,
	"RS6B" => 14,
	"RS6X" => 15,
	"RS7A" => 16,
	"RS7B" => 17,
	"RS7X" => 18,
	"RS8A" => 19,
	"RS8B" => 20,
	"RS8X" => 21,
	"RS9A" => 22,
	"RS9B" => 23,
	"RS9X" => 24,
	"RS10A" => 25,
	"RS10B" => 26,
	"RS10X" => 27,
	"RS11A" => 28,
	"RS11B" => 29,
	"RS11X" => 30,
	"RS12" => 31,
	"RS13" => 32,
	"RS14A" => 33,
	"RS14B" => 34,
	"RS14X" => 35,
	"RS15" => 36,
	"RS16A" => 37,
	"RS16B" => 38,
	"RS16X" => 39,
	"RS17A" => 40,
	"RS17B" => 41,
	"RS17X" => 42,
	"RS18A" => 43,
	"RS18B" => 44,
	"RS18X" => 45,
	"RS19A" => 46,
	"RS19B" => 47,
	"RS19X" => 48,
	"RS20" => 49,
	"RS21A" => 50,
	"RS21B" => 51,
	"RS21X" => 52,
	"RS22A" => 53,
	"RS22B" => 54,
	"RS22X" => 55,
	"RS23A" => 56,
	"RS23B" => 57,
	"RS23X" => 58,
	"RS24A" => 59,
	"RS24B" => 60,
	"RS24X" => 61,
	"RS25A" => 62,
	"RS25B" => 63,
	"RS25X" => 64,
	"RS26A" => 65,
	"RS26B" => 66,
	"RS26X" => 67,
	"RS27A" => 68,
	"RS27B" => 69,
	"RS27X" => 70,
	"RS28A" => 71,
	"RS28B" => 72,
	"RS28X" => 73,
	"RS29A" => 74,
	"RS29B" => 75,
	"RS29X" => 76,
	"RS30A" => 77,
	"RS30B" => 78,
	"RS30X" => 79,
	"RS31" => 80,
);

%pro2loc40Syeast = (
	"XXX" => 0,
	"S0A" => 1,
	"S0B" => 1,
	"S0X" => 1,
	"S1A" => 2,
	"S1B" => 2,
	"S1X" => 2,
	"S2" => 3,
	"S3" => 4,
	"S4A" => 5,
	"S4B" => 5,
	"S4X" => 5,
	"S5" => 6,
	"S6A" => 7,
	"S6B" => 7,
	"S6X" => 7,
	"S7A" => 8,
	"S7B" => 8,
	"S7X" => 8,
	"S8A" => 9,
	"S8B" => 9,
	"S8X" => 9,
	"S9A" => 10,
	"S9B" => 10,
	"S9X" => 10,
	"S10A" => 11,
	"S10B" => 11,
	"S10X" => 11,
	"S11A" => 12,
	"S11B" => 12,
	"S11X" => 12,
	"S12" => 13,
	"S13" => 14,
	"S14A" => 15,
	"S14B" => 15,
	"S14X" => 15,
	"S15" => 16,
	"S16A" => 17,
	"S16B" => 17,
	"S16X" => 17,
	"S17A" => 18,
	"S17B" => 18,
	"S17X" => 18,
	"S18A" => 19,
	"S18B" => 19,
	"S18X" => 19,
	"S19A" => 20,
	"S19B" => 20,
	"S19X" => 20,
	"S20" => 21,
	"S21A" => 22,
	"S21B" => 22,
	"S21X" => 22,
	"S22A" => 23,
	"S22B" => 23,
	"S22X" => 23,
	"S23A" => 24,
	"S23B" => 24,
	"S23X" => 24,
	"S24A" => 25,
	"S24B" => 25,
	"S24X" => 25,
	"S25A" => 26,
	"S25B" => 26,
	"S25X" => 26,
	"S26A" => 27,
	"S26B" => 27,
	"S26X" => 27,
	"S27A" => 28,
	"S27B" => 28,
	"S27X" => 28,
	"S28A" => 29,
	"S28B" => 29,
	"S28X" => 29,
	"S29A" => 30,
	"S29B" => 30,
	"S29X" => 30,
	"S30A" => 31,
	"S30B" => 31,
	"S30X" => 31,
	"S31" => 32,
);

%pro2loc60Syeast = (
	"XXX" => 0,
	"L1A" => 1,
	"L1B" => 1,
	"L1X" => 1,
	"L2A" => 2,
	"L2B" => 2,
	"L2X" => 2,
	"L3" => 3,
	"L4A" => 4,
	"L4B" => 4,
	"L4X" => 4,
	"L5" => 5,
	"L6A" => 6,
	"L6B" => 6,
	"L6X" => 6,
	"L7A" => 7,
	"L7B" => 7,
	"L7X" => 7,
	"L8A" => 8,
	"L8B" => 8,
	"L8X" => 8,
	"L9A" => 9,
	"L9B" => 9,
	"L9X" => 9,
	"L10" => 10,
	"L11A" => 11,
	"L11B" => 11,
	"L11X" => 11,
	"L12A" => 12,
	"L12B" => 12,
	"L12X" => 12,
	"L13A" => 13,
	"L13B" => 13,
	"L13X" => 13,
	"L14A" => 14,
	"L14B" => 14,
	"L14X" => 14,
	"L15A" => 15,
	"L15B" => 15,
	"L15X" => 15,
	"L16A" => 16,
	"L16B" => 16,
	"L16X" => 16,
	"L17A" => 17,
	"L17B" => 17,
	"L17X" => 17,
	"L18A" => 18,
	"L18B" => 18,
	"L18X" => 18,
	"L19A" => 19,
	"L19B" => 19,
	"L19X" => 19,
	"L20A" => 20,
	"L20B" => 20,
	"L20X" => 20,
	"L21A" => 21,
	"L21B" => 21,
	"L21X" => 21,
	"L22A" => 22,
	"L22B" => 22,
	"L22X" => 22,
	"L23A" => 23,
	"L23B" => 23,
	"L23X" => 23,
	"L24A" => 24,
	"L24B" => 24,
	"L24X" => 24,
	"L25" => 25,
	"L26A" => 26,
	"L26B" => 26,
	"L26X" => 26,
	"L27A" => 27,
	"L27B" => 27,
	"L27X" => 27,
	"L28" => 28,
	"L29" => 29,
	"L30" => 30,
	"L31A" => 31,
	"L31B" => 31,
	"L31X" => 31,
	"L32" => 32,
	"L33A" => 33,
	"L33B" => 33,
	"L33X" => 33,
	"L34A" => 34,
	"L34B" => 34,
	"L34X" => 34,
	"L35A" => 35,
	"L35B" => 35,
	"L35X" => 35,
	"L36A" => 36,
	"L36B" => 36,
	"L36X" => 36,
	"L37A" => 37,
	"L37B" => 37,
	"L37X" => 37,
	"L38" => 38,
	"L39" => 39,
	"L40A" => 40,
	"L40B" => 40,
	"L40X" => 40,
	"L41A" => 41,
	"L41B" => 41,
	"L41X" => 41,
	"L42A" => 42,
	"L42B" => 42,
	"L42X" => 42,
	"L43A" => 43,
	"L43B" => 43,
	"L43X" => 43,
	"P0" => 44,
	"P1A" => 45,
	"P1B" => 45,
	"P1X" => 45,
	"P2A" => 46,
	"P2B" => 46,
	"P2X" => 46,
);

%pro2loc80Syeast = (
	"XXX" => 0,
	"S0A" => 1,
	"S0B" => 1,
	"S0X" => 1,
	"S1A" => 2,
	"S1B" => 2,
	"S1X" => 2,
	"S2" => 3,
	"S3" => 4,
	"S4A" => 5,
	"S4B" => 5,
	"S4X" => 5,
	"S5" => 6,
	"S6A" => 7,
	"S6B" => 7,
	"S6X" => 7,
	"S7A" => 8,
	"S7B" => 8,
	"S7X" => 8,
	"S8A" => 9,
	"S8B" => 9,
	"S8X" => 9,
	"S9A" => 10,
	"S9B" => 10,
	"S9X" => 10,
	"S10A" => 11,
	"S10B" => 11,
	"S10X" => 11,
	"S11A" => 12,
	"S11B" => 12,
	"S11X" => 12,
	"S12" => 13,
	"S13" => 14,
	"S14A" => 15,
	"S14B" => 15,
	"S14X" => 15,
	"S15" => 16,
	"S16A" => 17,
	"S16B" => 17,
	"S16X" => 17,
	"S17A" => 18,
	"S17B" => 18,
	"S17X" => 18,
	"S18A" => 19,
	"S18B" => 19,
	"S18X" => 19,
	"S19A" => 20,
	"S19B" => 20,
	"S19X" => 20,
	"S20" => 21,
	"S21A" => 22,
	"S21B" => 22,
	"S21X" => 22,
	"S22A" => 23,
	"S22B" => 23,
	"S22X" => 23,
	"S23A" => 24,
	"S23B" => 24,
	"S23X" => 24,
	"S24A" => 25,
	"S24B" => 25,
	"S24X" => 25,
	"S25A" => 26,
	"S25B" => 26,
	"S25X" => 26,
	"S26A" => 27,
	"S26B" => 27,
	"S26X" => 27,
	"S27A" => 28,
	"S27B" => 28,
	"S27X" => 28,
	"S28A" => 29,
	"S28B" => 29,
	"S28X" => 29,
	"S29A" => 30,
	"S29B" => 30,
	"S29X" => 30,
	"S30A" => 31,
	"S30B" => 31,
	"S30X" => 31,
	"S31" => 32,

	"L1A" => 33,
	"L1B" => 33,
	"L1X" => 33,
	"L2A" => 34,
	"L2B" => 34,
	"L2X" => 34,
	"L3" => 35,
	"L4A" => 36,
	"L4B" => 36,
	"L4X" => 36,
	"L5" => 37,
	"L6A" => 38,
	"L6B" => 38,
	"L6X" => 38,
	"L7A" => 39,
	"L7B" => 39,
	"L7X" => 39,
	"L8A" => 40,
	"L8B" => 40,
	"L8X" => 40,
	"L9A" => 41,
	"L9B" => 41,
	"L9X" => 41,
	"L10" => 42,
	"L11A" => 43,
	"L11B" => 43,
	"L11X" => 43,
	"L12A" => 44,
	"L12B" => 44,
	"L12X" => 44,
	"L13A" => 45,
	"L13B" => 45,
	"L13X" => 45,
	"L14A" => 46,
	"L14B" => 46,
	"L14X" => 46,
	"L15A" => 47,
	"L15B" => 47,
	"L15X" => 47,
	"L16A" => 48,
	"L16B" => 48,
	"L16X" => 48,
	"L17A" => 49,
	"L17B" => 49,
	"L17X" => 49,
	"L18A" => 50,
	"L18B" => 50,
	"L18X" => 50,
	"L19A" => 51,
	"L19B" => 51,
	"L19X" => 51,
	"L20A" => 52,
	"L20B" => 52,
	"L20X" => 52,
	"L21A" => 53,
	"L21B" => 53,
	"L21X" => 53,
	"L22A" => 54,
	"L22B" => 54,
	"L22X" => 54,
	"L23A" => 55,
	"L23B" => 55,
	"L23X" => 55,
	"L24A" => 56,
	"L24B" => 56,
	"L24X" => 56,
	"L25" => 57,
	"L26A" => 58,
	"L26B" => 58,
	"L26X" => 58,
	"L27A" => 59,
	"L27B" => 59,
	"L27X" => 59,
	"L28" => 60,
	"L29" => 61,
	"L30" => 62,
	"L31A" => 63,
	"L31B" => 63,
	"L31X" => 63,
	"L32" => 64,
	"L33A" => 65,
	"L33B" => 65,
	"L33X" => 65,
	"L34A" => 66,
	"L34B" => 66,
	"L34X" => 66,
	"L35A" => 67,
	"L35B" => 67,
	"L35X" => 67,
	"L36A" => 68,
	"L36B" => 68,
	"L36X" => 68,
	"L37A" => 69,
	"L37B" => 69,
	"L37X" => 69,
	"L38" => 70,
	"L39" => 71,
	"L40A" => 72,
	"L40B" => 72,
	"L40X" => 72,
	"L41A" => 73,
	"L41B" => 73,
	"L41X" => 73,
	"L42A" => 74,
	"L42B" => 74,
	"L42X" => 74,
	"L43A" => 75,
	"L43B" => 75,
	"L43X" => 75,
	"P0" => 76,
	"P1A" => 77,
	"P1B" => 77,
	"P1X" => 77,
	"P2A" => 78,
	"P2B" => 78,
	"P2X" => 78,
);

%pro2locRNA = (
	"XXX" => 0,
	"16S" => 1,
	"23S" => 2,
	"5s" => 3,
);

1;
