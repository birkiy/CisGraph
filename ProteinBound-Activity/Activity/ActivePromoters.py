


Mat = pd.read_csv("Concat.tab", sep="\t")

import math
# mRNA for TSS +

highP = Mat[Mat["deepTools_group"] == "tsP.SecondPart.bed"]
highM = Mat[Mat["deepTools_group"] == "tsM.FirstPart.bed"]

highPactivity = highP["groseq_dht.plus"] / ((highP["groseq_dht.minus"] + highP["groseq_dht.minus.1"] + 1 )/2) + 1


highMactivity = highM["groseq_dht.minus"] / ((highM["groseq_dht.plus"] + highM["groseq_dht.plus.1"] + 1 )/2) + 1

highP["Activity"] = np.log(highPactivity)
highM["Activity"] = np.log(highMactivity)


tsP = highP[["#chrom", "start", "end", "name", "Activity"]]
tsP["start"] -= 500

tsP.to_csv("tsPActivity.bed", sep="\t", index=False, header=False)

tsM = highM[["#chrom", "start", "end", "name", "Activity"]]
tsM["end"] += 500

tsM.to_csv("tsMActivity.bed", sep="\t", index=False, header=False)





Mat = pd.read_csv("ConcatARBS.tab", sep="\t")

import math
# mRNA for TSS +

conP = Mat[Mat["deepTools_group"] == "con.ARBS.P.bed"] # Plus Second
conP["Strand"] = "+"
conP["nodeClass"] = "con"
conP["start"] -= 350
conM = Mat[Mat["deepTools_group"] == "con.ARBS.M.bed"] # Minus First
conM["Strand"] = "-"
conM["nodeClass"] = "con"
conM["end"] += 350

indP = Mat[Mat["deepTools_group"] == "ind.ARBS.P.bed"] # Plus Second
indP["Strand"] = "+"
indP["nodeClass"] = "ind"
indP["start"] -= 350
indM = Mat[Mat["deepTools_group"] == "ind.ARBS.M.bed"] # Minus First
indM["Strand"] = "-"
indM["nodeClass"] = "ind"
indM["end"] += 350

nonP = Mat[Mat["deepTools_group"] == "non.ARBS.P.bed"] # Plus Second
nonP["Strand"] = "+"
nonP["nodeClass"] = "non"
nonP["start"] -= 350
nonM = Mat[Mat["deepTools_group"] == "non.ARBS.M.bed"] # Minus First
nonM["Strand"] = "-"
nonM["nodeClass"] = "non"
nonM["end"] += 350




highP = conP
highM = conM

highPactivity = highP["groseq_dht.plus"] / ((highP["groseq_dht.minus"] + highP["groseq_dht.minus.1"] + 1 )/2) + 1

highMactivity = highM["groseq_dht.minus"] / ((highM["groseq_dht.plus"] + highM["groseq_dht.plus.1"] + 1 )/2) + 1

highP["Activity"] = np.log(highPactivity)
highM["Activity"] = np.log(highMactivity)

conP = highP
conM = highM


highP = indP
highM = indM

highPactivity = highP["groseq_dht.plus"] / ((highP["groseq_dht.minus"] + highP["groseq_dht.minus.1"] + 1 )/2) + 1

highMactivity = highM["groseq_dht.minus"] / ((highM["groseq_dht.plus"] + highM["groseq_dht.plus.1"] + 1 )/2) + 1

highP["Activity"] = np.log(highPactivity)
highM["Activity"] = np.log(highMactivity)

indP = highP
indM = highM




highP = nonP
highM = nonM

highPactivity = highP["groseq_dht.plus"] / ((highP["groseq_dht.minus"] + highP["groseq_dht.minus.1"] + 1 )/2) + 1

highMactivity = highM["groseq_dht.minus"] / ((highM["groseq_dht.plus"] + highM["groseq_dht.plus.1"] + 1 )/2) + 1

highP["Activity"] = np.log(highPactivity)
highM["Activity"] = np.log(highMactivity)

nonP = highP
nonM = highM



arbsArctivity = pd.concat([conP, conM, indP, indM, nonP, nonM])

arbsArctivity.to_csv("arbsActivity.bed", sep="\t", index=False, header=False)
