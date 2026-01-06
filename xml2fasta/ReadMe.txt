unirefxml2fasta.pl
------------------------------------------
A perl script for quick conversion from UniRef xml files to fasta files. The result will be gzipped automatically.
URL: https://proteininformationresource.org/download/uniref/xml2fasta/unirefxml2fasta.pl

Usage:
perl unirefxml2fasta.pl uniref90.xml.gz new_uniref90.fasta.gz > convert.log

The output file name is automatically generated if not present.

*Note:
There has been a format change for UniRef xml files since 2019_11 (same as other UniProt xml files) where sequence
is presented differently. This script relies on the release version (e.g. releaseDate="2020-04-22" version="2020_02")
in the xml to automatically handle the version difference. If you run the script on partial xml please make sure it
carries said information.

-------------------------------------------
License: GPLv3, https://www.gnu.org/licenses/gpl%2D3.0.html
