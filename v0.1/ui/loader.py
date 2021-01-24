import json
def loadData(req):
    f = open('promod/data/analysis.json','r')
    analysis_data = json.loads(f.read())
    if len(req) == 1:
        found = False
        print('Files which are available to run in analysis.json:')
        for file in range(0,len(analysis_data)):
            filename = analysis_data[file]['file_name']
            print(f'\t{str(file)} - {filename}')
        while not found:
            isInt = False
            isFile = False
            n = input("\nPlease enter the number of the saved analysis you'd like to run\t")
            try:
                n = int(n)
                isInt = True
            except:
                print('\nYou must enter an integer')
            if isInt:
                try:
                    data = analysis_data[n]
                    isFile = True
                except:
                    print('\nNo file associated with that number')
            if isInt and isFile:
                found = True
                print(f'\nRunning analysis {n}, please wait...\n')
        return data
    elif len(req) == 2:
        if (req[1] == 'stdin'):
            data = json.loads(sys.stdin.read())
        else:
            isInt = False
            try:
                n = int(req[1])
                isInt = True
            except:
                print('\nTo load a file you must use the integer associated with that file\n')
                return 
            if isInt:
                try:
                    data = analysis_data[n]
                    print(f'\nRunning analysis {n}, please wait...\n')
                except:
                    print('\nNo file associated with that number\n')
                    return
                return data
    else:
        print('\nWrong number of arguments\n')
        return