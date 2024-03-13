#Codigo para optener solo los datos de T2M(temperatura) y selecionar el dato maximo por dia 
#!/bin/bash
vars="daymax"
for var in $vars
do
 for f in `ls -d *`
  do
  cdo -b f64 addc,-273.15 -selname,t2m -${var} ${f:0} max/${f:0:-3}_${var}_t2m_c.nc
  done
done 

#Codigo para optener solo los datos de T2M(temperatura) y selecionar el dato minimo por dia 
#!/bin/bash
vars="daymin"
for var in $vars
do
 for f in `ls -d *`
  do
  cdo -b f64 addc,-273.15 -selname,t2m -${var} ${f:0} min/${f:0:-3}_${var}_t2m_c.nc
  done
done 