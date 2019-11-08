from collections import defaultdict
from shapely.strtree import STRtree
import shapely
import shapely.wkt
from shapely.ops import transform
import shapely.geometry
import os
import sys

#INPUT_FILE_NAME="test.tsv"
INPUT_FILE_NAME="ops_ready_boundaries.tsv"
ERROR_FILE_NAME="errors.txt"
OUTPUT_FILE_NAME="overlaps.txt"

SIM_THRESH=0.05

fp = open(INPUT_FILE_NAME, 'r')
#File format: a set of boundaries to check for overlap - if check_overlap is True then we will check that boundary for overlap against all others in the list. 
#Assumes the first line is a header. 
#POLY_WKT FIELD_ID CHECK_OVERLAP?
num_fields=3
poly_ind=0
field_ind=1
check_ind=2

# Using query (but not in this script
query=f"""select A.field_wkt, A.dca_field_id fid, case when B.status='OPS_READY' then TRUE else FALSE end as is_ops_ready
from 
PARTHENON.AGRONOMY.DCA_GROWER_FARM_FIELD A 
left join (select * from PARTHENON.DCA_DB.DCA_FIELD_BOUNDARYS where dca_boundary_id in (select active_boundary_id from  PARTHENON.DCA_DB.DCA_FIELDS)
and is_active='TRUE') B
on A.dca_field_id=B.dca_field_id"""


errors = open(ERROR_FILE_NAME, 'w')

#To be populated in this loop. 
geo_index={} #A mapping from WKT string to the input line so that we can go from geo back to field. We use the WKT produced by the shapely library. It actually maps to a list because we could have multiple fields with exactly the same WKT
all_geos=[] #A list of all of the geometries in the file, this will be what we compare against 
check_geos=[] #A list of the geometries in the file we want to check, each entry  [input_field0, input_field1, ... input_field_n, POLY obj]
check_field_ids=set() #Keep a set of all of the field IDs that we check to avoid double counting

header=True
# Read the input files to get the boundaries
for line in fp:
    #Skip the header line
    if header:
        header=False
        continue
    
    words = line.strip().split('\t')
    words[check_ind] = (words[check_ind] == 'TRUE') #translate to bool

    if len(words) != num_fields:
        print(f"Input file not correctly formatted. Should have {num_fields} fields per line")
        print(line)
        exit(1)

    try: 
        obj = shapely.wkt.loads(words[poly_ind])
        if not obj.is_valid:
            type = obj.geom_type
            errors.write(f"Invalid polygon of type {type} field_id {words[field_ind]} check? {words[check_ind]} WKT:\n{line}\n")
            continue
    except:
        errors.write(f"Error parsing geometry field_id {words[field_ind]} check? {words[check_ind]} WKT:\n{line}\n")
    else:

        #Populate the geometry lists. 
        all_geos.append(obj)
        if words[check_ind] or True:
            words.append(obj)
            check_geos.append(words)
        wkt_b = shapely.wkt.dumps(obj)

        if wkt_b not in geo_index:
            geo_index[wkt_b] = [words]
        else:
            geo_index[wkt_b].append(words)

fp.close()

print(f"Done processing boundaries. Total {len(all_geos)} To check {len(check_geos)}")

#Create Tree to store all the geos to compare against
tree = STRtree(all_geos)

outf = open(OUTPUT_FILE_NAME, 'w')
overlap_count=0

#Iterate through the geometries we want to check overlap
fieldpairset = set()

outf.write(f"Similarity;FieldID_A;FieldID_B;WKT_A;WKT_B\n")
for item1 in check_geos:
    geo1 = item1[num_fields] #This is the last field, we appended the polygon object. 
    overlaps = tree.query(geo1)
    for geo2 in overlaps:
        try:
            intersection_size = geo1.intersection(geo2).area
            union_size = geo1.union(geo2).area
            sim = intersection_size/union_size
            if sim > SIM_THRESH:
                field1 = item1[field_ind]
                #find field2
                ind = shapely.wkt.dumps(geo2)
                if ind in geo_index:
                    field2_list = geo_index[ind]
                else:
                    print("We should never get here. Every geo should have a corresponding field")
                for item2 in field2_list:
                    field2=item2[field_ind]
                    if field1 == field2:
                        continue
                    #Avoid double counting
                    pairstr = f"{min(field1,field2)} {max(field1,field2)}"
                    if pairstr in fieldpairset:
                        continue
                    else:
                        fieldpairset.add(pairstr)

                    outf.write(f"{sim};{field1};{field2};{item1[poly_ind]};{item2[poly_ind]}\n")
                    overlap_count+=1

        except shapely.errors.TopologicalError:
            errors.write(f"Unexpected error: {sys.exc_info()[0]}")
            errors.write(shapely.wkt.dumps(geo1))
            errors.write(shapely.wkt.dumps(geo2))

outf.close()
errors.close()
print(f"Overlap Count {overlap_count}")
