
seq = input("Insert sequence: ")

seq = seq.upper()

puri = seq.count("A") + seq.count("G")
pyri = seq.count("C") + seq.count("T")

print("Purines: ", puri)
print("Pyrimidines: ", pyri)