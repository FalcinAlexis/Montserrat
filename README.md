# Repository for sharing codes for getting Montserrat Seisan database into ObsPy-CSV files

# Step 1: setup a .bashrc_montserrat file in your ~ (HOME) directory
# Example of file contents:
export PYTHONPATH="/Users/thompsong/src/Alexis_Montserrat_codes"
export PATH=$PATH:/Users/thompsong/src/Alexis_Montserrat_codes
export SEISAN_DATA="/Users/thompsong/Desktop/seismo"

# Step 2: source this from your .bashrc by doing
echo "source ~/.bashrc_montserrat" >> ~/.bashrc

# Step 3: start a new bash terminal

# Step 4: run seisandb2csv.py
cd $SEISAN_DATA
seisandb2csv.py

# Step 5: A CSV file like MVOE_catalog200501.csv should now have been created for each month that S-files are in your Seisan database, e.g. under $SEISAN_DATA/REA. You could now use standard Linux commands, e.g.

# (a) grep out 1 channel only
grep MBWH.Z.BH MVOE_catalog*.csv > MBWH.BHZ.csv

# (b) grep out 1 channel and count duration 
grep MBWH.Z.BH MVOE_catalog*.csv | awk -F',' '{print $4}' 

# (c) count the number of each volcanic subclass
awk -F',' '{print $3}' MVOE_*.csv | sort | uniq -c

# (d) count the number of seconds total for each different channel
awk -F',' '{a[$9]+=$4}END{for(i in a) print i,a[i]}' MV*csv | sort

# Step 6: CSV files can directly be imported in Python dataframes and manipulated there.

# Step 7: Alexis, this would be a great place to also put any *.pynb and *.json and other *.py codes related to running malfante/AAA codes on the Montserrat data

Glenn Thompson 2019/10/16, at IPGP
