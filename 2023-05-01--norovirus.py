group_I = 0
group_II = 0
group_IV = 0
group_IX = 0
group_I_II = 0
group_II_IV = 0
group_I_II_IX = 0
group_any = 0


with open("prevalence-data/cdc-nors-outbreak-data.tsv") as inf:
    cols = None
    for line in inf:
        row = line.strip().split("\t")
        if not cols:
            cols = row
            continue
        etiologies = set(row[cols.index("Etiology")].split("; "))
        etiologies = set(x.replace(" unknown", "").replace(" other", "")
                         for x in etiologies
                         if "Norovirus" in x)
        if len(etiologies) > 1 and "Norovirus" in etiologies:
            etiologies.remove("Norovirus")
        
        genotype = row[cols.index("Serotype or Genotype")]

        if not any("Norovirus" in etiology for etiology in etiologies):
            continue

        group_any += 1
        if len(etiologies) == 1:
            if list(etiologies) == ['Norovirus Genogroup I']:
                group_I += 1
            elif list(etiologies) == ['Norovirus Genogroup II']:
                group_II += 1
            elif list(etiologies) == ['Norovirus Genogroup IV']:
                group_IV += 1
            elif list(etiologies) == ['Norovirus Genogroup IX']:
                group_IX += 1
            elif list(etiologies) == ['Norovirus']:
                pass
            else:
                raise Exception(etiologies)
        elif len(etiologies) == 2:
            if etiologies == set(('Norovirus Genogroup I',
                                  'Norovirus Genogroup II')):
                group_I_II += 1
            elif etiologies == set(('Norovirus Genogroup II',
                                    'Norovirus Genogroup IV')):
                group_II_IV += 1
            else:
                raise Exception(etiologies)
        elif len(etiologies) == 3:
            if etiologies == set(('Norovirus Genogroup I',
                                  'Norovirus Genogroup II',
                                  'Norovirus Genogroup IX')):
                group_I_II_IX += 1
            else:
                raise Exception(etiologies)
        else:
            raise Exception(etiologies)

        #print ("%s\t%s" % (";".join(etiologies), genotype))
        #print (etiologies)

print("I\t%s" % (group_I))
print("II\t%s" % (group_II))
print("IV\t%s" % (group_IV))
print("IX\t%s" % (group_IX))
print("I+II\t%s" % (group_I_II))
print("II+IV\t%s" % (group_II_IV))
print("I+II+IX\t%s" % (group_I_II_IX))
print(group_any)

