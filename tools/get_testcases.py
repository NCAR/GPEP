import argparse, os, sys, pathlib, requests

# parse arguments
parser = argparse.ArgumentParser(description='Parse input arguments')
parser.add_argument('-n', '--number', type=str, help='The number of test cases that will be downloaded')
parser.add_argument('-o', '--outpath', type=str, help='The output path of the test case')
args = parser.parse_args()

if args.number is not None:
    print(f'The value of -n option is: {args.number}')
else:
    print('All test cases will be downloaded')
    
if  args.outpath is not None:
    outpath = args.outpath
    print(f'The output path is :{args.outpath}')
else:
    path0 = str(pathlib.Path(os.path.abspath(sys.argv[0])).parents[2])
    outpath = path0 + '/GPEP_test_cases'
    print(f'Use default output path: {outpath}')
    
os.makedirs(outpath, exist_ok=True)

# download test cases
urls = [ ['https://zenodo.org/record/8222852/files/cali2017.tar.gz',  'cali2017.tar.gz'],
         ['https://zenodo.org/record/8222852/files/UpCO_SWE.tar.gz',  'UpCO_SWE.tar.gz'],
       ]

if args.number is None:
    numbers = [0]
else:
    numbers = [int(i) for i in args.number.split(',')]

for i in numbers:
    urli = urls[i][0]
    outfile = outpath + '/' + urls[i][1]
    print(f'Download {urli} to {outfile}')
    
    # Download
    r = requests.get(urli, allow_redirects=True)
    open(outfile, 'wb').write(r.content)
    
    # uncompress
    print('Uncompress files')
    cwd = os.getcwd()
    os.chdir(outpath)
    if urls[i][1].endswith('tar') or urls[i][1].endswith('tar.gz'): 
        os.system(f'tar -xf {urls[i][1]}')
    elif urls[i][1].endswith('zip'):
        os.system(f'unzip {urls[i][1]}')
    else:
        print('Unsupported format')
    os.chdir(cwd)
    