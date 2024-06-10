
def config():
    print("Do you want to use your most recent project name, genome, design matrix, and reads folders:(e.g y/n)")
    print("\033[91m"+"If this is your first time, type"+"\033[94m"+" n"+"\x1b[0m")
    param = input()
    if param == 'n':
        print("========================================")
        print("Provide any unique project name:")
        print("\033[91m"+"If you want to use our example dataset, type"+"\033[94m"+" example_dataset"+"\x1b[0m")
        project_name = input()
        conf = open("../scripts_DoNotTouch/config.py", "w")
        conf.write("project_name="+"'"+project_name+"'"+"\n")
        conf.write("parameters_exist="+"'"+param+"'"+"\n")
        conf.close()
    else:
        with open('../scripts_DoNotTouch/config.py', 'r') as file:
            data = file.readlines()
        data[1] = "parameters_exist="+"'"+param+"'"+"\n"
        with open('../scripts_DoNotTouch/config.py', 'w') as file:
            file.writelines( data )
            
