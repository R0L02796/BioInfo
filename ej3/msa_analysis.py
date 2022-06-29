
print("Gorilla")
gorilla_homo = open("outputs/clustal/gorilla_homo", "r")
print("Secuencias con similitudes: ", len(gorilla_homo.readlines()))
gorilla_homo = open("outputs/clustal/gorilla_homo", "r").read()
print("Matches totales: ", gorilla_homo.count("*"))
print("Matches con score >= 0.5 (match mayor al 50%) : ", gorilla_homo.count(":"))
print("Matches con score < 0.5 (match menor al 50%) : ", gorilla_homo.count("."))

print("\n")

print("Pan homo (Chimpance)")
pan_homo = open("outputs/clustal/pan_homo", "r")
print("Secuencias con similitudes: ", len(pan_homo.readlines()))
pan_homo = open("outputs/clustal/pan_homo", "r").read()
print("Matches totales: ", pan_homo.count("*"))
print("Matches con score >= 0.5 (match mayor al 50%) : ", pan_homo.count(":"))
print("Matches con score < 0.5 (match menor al 50%) : ", pan_homo.count("."))
print("\n")

print("Vulpes (Zorro artico)")
vulpes_homo = open("outputs/clustal/vulpes_homo", "r")
print("Secuencias con similitudes: ", len(vulpes_homo.readlines()))
vulpes_homo = open("outputs/clustal/vulpes_homo", "r").read()
print("Matches totales: ", vulpes_homo.count("*"))
print("Matches con score >= 0.5 (match mayor al 50%) : ", vulpes_homo.count(":"))
print("Matches con score < 0.5 (match menor al 50%) : ", vulpes_homo.count("."))

print("\n")
