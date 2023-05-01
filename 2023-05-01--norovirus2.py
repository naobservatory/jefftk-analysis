data = {} # year -> [all, has etiology breakdown]

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

        group = None
        if len(etiologies) == 1:
            if list(etiologies) == ['Norovirus Genogroup I']:
                group = "I"
            elif list(etiologies) == ['Norovirus Genogroup II']:
                group = "II"
            elif list(etiologies) == ['Norovirus Genogroup IV']:
                group = "IV"
            elif list(etiologies) == ['Norovirus Genogroup IX']:
                group = "IX"
            elif list(etiologies) == ['Norovirus']:
                pass
            else:
                raise Exception(etiologies)
        elif len(etiologies) == 2:
            if etiologies == set(('Norovirus Genogroup I',
                                  'Norovirus Genogroup II')):
                group = "I+II"
            elif etiologies == set(('Norovirus Genogroup II',
                                    'Norovirus Genogroup IV')):
                group = "II+IV"
            else:
                raise Exception(etiologies)
        elif len(etiologies) == 3:
            if etiologies == set(('Norovirus Genogroup I',
                                  'Norovirus Genogroup II',
                                  'Norovirus Genogroup IX')):
                group = "I+II+IX"
            else:
                raise Exception(etiologies)
        else:
            raise Exception(etiologies)

        year = row[cols.index("Year")]

        if year not in data:
            data[year] = [0,0]

        data[year][0] += 1
        if group:
            data[year][1] += 1

for year in sorted(data):
    print("%s\t%s\t%s" % (
        year, *data[year]))
            
        

        

