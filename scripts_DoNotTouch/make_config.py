
def config():
    
    print("Do you want to use your most recent project name, genome, design matrix, and reads folders:(e.g y/n)")
    print("\033[91m"+"If this is your first time, type"+"\033[94m"+" n"+"\x1b[0m")
    param = input()
    if param == 'n':
        print("========================================")
        print("Provide any unique project name:")
        print("\033[91m"+"If you want to use our example dataset, type"+"\033[94m"+" example_dataset"+"\x1b[0m")
        project_name = input()

        print("========================================")
        print("Provide a path/location where to store your analysis results:(Sometimes your home has limited space)")
        print("\033[91m"+"If not sure, leave it blank and press enter. Results will be stored in the same location"+"\033[94m")
        res_dir = input()
        if len(res_dir) == 0:
            res_dir = "../../csl_results/"
        else:
            res_dir = res_dir+"/"
        
        conf = open("../scripts_DoNotTouch/config.py", "w")
        conf.write("project_name="+"'"+project_name+"'"+"\n")
        conf.write("parameters_exist="+"'"+param+"'"+"\n")
        conf.write("results_directory="+"'"+res_dir+"'"+"\n")
        conf.close()
    else:
        with open('../scripts_DoNotTouch/config.py', 'r') as file:
            data = file.readlines()
        data[1] = "parameters_exist="+"'"+param+"'"+"\n"
        with open('../scripts_DoNotTouch/config.py', 'w') as file:
            file.writelines( data )
            
