enzyme_file_name = input("Input enzyme file name: ")
protein_ligand_file_name = input("Input protein ligand file name: ")
complex_name = input("Input complex name: ")

for i in range(1, int(input("How many ligand files do you have?: "))  + 1):
	with open(f"{complex_name}_{i}.pdb", "w+") as f:
		with open(f"{enzyme_file_name}.pdb", "r") as e:
			f.write("".join(e.readlines()))
		with open(f"{protein_ligand_file_name}_{i}.pdb", "r") as l:
			for line in l:
				if "ATOM" in line:
					line_list = list(line)
					line_list[21] = "X"
					line = "".join(line_list)
					f.write(line)
				f.write(line)



