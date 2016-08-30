#! /usr/bin/python
#
# join
#   Joing pages from a a collection of PDF files into a single PDF file.
#
#   join [--output <file>] [--append] [--shuffle] [--preview] [--verbose]"
#
#   Parameter:
#
#   --shuffle
#	Take a page from each PDF input file in turn before taking another from each file.
#	If this option is not specified then all of the pages from a PDF file are appended
#	to the output PDF file before the next input PDF file is processed.
#
#   --verbose
#   Write information about the doings of this tool to stderr.
#
import sys
import os
import getopt
import tempfile
import shutil
from CoreGraphics import *

def createPDFDocumentWithPath(path):
	print path
	return CGPDFDocumentCreateWithProvider(CGDataProviderCreateWithFilename(path))
	
def shufflePages(writeContext, inputFiles, verbose):
	
	maxPages = 0

	print inputFiles
	
	# create PDFDocuments for all of the files.
	docs = map(createPDFDocumentWithPath, inputFiles)

	print docs
	
	# find the maximum number of pages.
	for doc in docs:
		if doc.getNumberOfPages() > maxPages:
			maxPages = doc.getNumberOfPages()
			print maxPages
			
	# Shuffle the pages
	for pageNum in xrange(1, maxPages + 1):
		print "pageNum %d" % (pageNum)
		for doc in docs:
			page = doc.getPage(pageNum)
			if page != None:
				mediaBox = doc.getMediaBox(pageNum)
				writeContext.beginPage(mediaBox)
				writeContext.drawPDFDocument(mediaBox, doc, pageNum)
				writeContext.endPage()
				
def append(writeContext, inputFiles, verbose):
	"""Write the pages from supplied PDF files into the provided context"""
	for file in inputFiles:
		inputFile = CGPDFDocumentCreateWithProvider(CGDataProviderCreateWithFilename(file))
		for pageNum in xrange(1, inputFile.getNumberOfPages() + 1) :
			mediaBox = inputFile.getMediaBox(pageNum)
			writeContext.beginPage(mediaBox)
			writeContext.drawPDFDocument(mediaBox, inputFile, pageNum)
			writeContext.endPage()
			if verbose:
				print "Copied page %d from %s" % (pageNum, file)

def main(argv):

	# The PDF context we will draw into to create a new PDF
	writeContext = None

	# If True then generate more verbose information
	verbose = False
	source = None
	shuffle = False
	
	# Parse the command line options
	try:
		options, args = getopt.getopt(argv, "o:sv", ["output=", "shuffle", "verbose"])

	except getopt.GetoptError:
		usage()
		sys.exit(2)

	for option, arg in options:

		if option in ("-o", "--output") :
			if verbose:
				print "Setting %s as the destination." % (arg)
			writeFilename = arg
			pageRect = CGRectMake (0, 0, 612, 792)
			writeContext = CGPDFContextCreateWithFilename(arg, pageRect)

		elif option in ("-s", "--shuffle") :
			if verbose :
				print "Shuffle pages to the output file."
			shuffle = True

		elif option in ("-v", "--verbose") :
			print "Verbose mode enabled."
			verbose = True

		else :
			print "Unknown argument: %s" % (option)

	if shuffle:
		shufflePages(writeContext, args, verbose)
	else:
		append(writeContext, args, verbose)
    
def usage():
	print "Usage: join [--output <file>] [--append] [--shuffle] [--preview] [--verbose]"

if __name__ == "__main__":
	main(sys.argv[1:])

    