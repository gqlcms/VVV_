# 18 code
# /BulkGravToWW_narrow_M-1000_TuneCP5_PSWeights_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM
import os
num_old_string = ["template_1","template_2"]
template_file = "template.py"
submit = ""
check = ""
resubmit =""
os.system("mkdir -p 17crabfolder")
for i in (open("17samplelist.txt", 'r').readlines()):
     dataset = i.replace("\n","")
     filename = i.split("/")
     if len(filename[1]) >0 :       
        if "WWZ_4F_" in filename[1]:
            newfilatag = "17_"+filename[1].split("_TuneCP5")[0]
	    newfilename = filename[1].split("_TuneCP5")[0]+".py"
            submit =  submit+ "crab submit -c ./17code/17crabfolder/" + newfilename+"\n"
            check = check+ "crab status -d ./crab_"+newfilatag+"\n"
            resubmit = resubmit+ "crab resubmit ./crab_"+newfilatag+"\n"
            print newfilename            
            file_data = ""
            with open( template_file , 'r') as f:
                for line in f:
                    if "template_1" in line :
                       line = line.replace( "template_1" , newfilatag)
                    if "template_2" in line :
                       line = line.replace( "template_2" , dataset )
                    file_data += line
            with open("17crabfolder/"+newfilename,"w") as f:
                f.write(file_data)
with open("../submit_17.sh","w") as f:
     f.write(submit)
with open("../check_17.sh","w") as f:
     f.write(check)
with open("../resubmit_17.sh","w") as f:
     f.write(resubmit)
os.system("chmod 755 ../submit_17.sh;chmod 755 ../check_17.sh;chmod 755 ../resubmit_17.sh")

