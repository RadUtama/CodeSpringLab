
def config():
    print("Provide any unique project name:")
    print("\033[91m"+"If you want to use our example dataset, type"+"\033[94m"+" example_dataset"+"\x1b[0m")
    project_name = input()
    conf = open("../scripts_DoNotTouch/config.py", "w")
    conf.write("project_name="+"'"+project_name+"'")
    conf.close()