@echo off
echo Downloading BLAST+ installer...
powershell -Command "(New-Object Net.WebClient).DownloadFile('ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.7.1/ncbi-blast-2.7.1+-win64.exe', 'ncbi-blast-2.7.1+-win64.exe')"
echo Running BLAST+ installer...
start /wait ncbi-blast-2.7.1+-win64.exe /S
echo Testing BLAST+ installation...
"C:\Program Files\NCBI\blast-2.7.1+\bin\blastp.exe" -h
pause