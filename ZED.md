---
title: "Projekt ZED"
author: "Piotr Skoczek"
date: "2 grudnia 2018"
output: 
  html_document: 
    keep_md: yes
    number_sections: yes
    toc: yes
    toc_float: yes
---



#Wprowadzenie

Niniejszy raport przedstawia analizê danych w jêzyku R. Zestaw danych, pochodz¹cy z Protein Data Bank, dotyczy ligandów. Po wstêpnym przetworzeniu danych pokazano rozk³ad liczby atomów i elektronów, ligandy o najwiêkszej niezgodnoœci liczby atomów i elektronów, czy korelacjê miêdzy niektórymi zmiennymi. Raport zakoñczony jest stworzeniem klasyfikatora, który na podstawie pozosta³ych parametrów próbuje przewidzieæ nazwê ligandu.

#£adowanie bibliotek


```r
library(data.table)
library(dplyr)
library(DT)
library(tidyr)
library(ggplot2)
library(plotly)
library(knitr)
library(reshape2)
library(kableExtra)
library(caret)
library(party)
opts_chunk$set(message = FALSE, warning = FALSE)
```

#Zapewnienie powtarzalnoœci


```r
set.seed(123)
```

#Wczytywanie danych z pliku


```r
data <- fread("all_summary.csv", nrows = 300000)
```


#Usuwanie niektórych wierszy


```r
to_remove <- c("UNK", "UNX", "UNL", "DUM", "N", "BLOB", "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "MSE", "PHE", "PRO", "SEC", "SER", "THR", "TRP", "TYR", "VAL", "DA", "DG", "DT", "DC", "DU", "A", "G", "T", "C", "U", "HOH", "H20", "WAT")
data_cleared <- data %>% select(-(blob_coverage:pdb_code), -(res_id:local_res_atom_count), -(local_res_atom_non_h_occupancy_sum), -(local_res_atom_non_h_electron_occupancy_sum:local_res_atom_S_count), -(dict_atom_C_count:skeleton_periphery), -(local_max_over_std), -(local_cut_by_mainchain_volume:local_near_cut_count_N), -(fo_col:resolution_max_limit), -(part_step_FoFc_std_min:part_step_FoFc_std_step)) %>%
  filter(res_name != to_remove)
```

Powy¿szy kod usuwa wiersze z podanymi wartoœciami res_name oraz kolumny, które s¹ nieopisane, niewykorzystywane do klasyfikacji lub zawieraj¹ce du¿o brakuj¹cych wartoœci.

#Usuwanie wierszy z brakuj¹cymi wartoœciami


```r
data_final <- drop_na(data_cleared)
```

#Podsumowanie

Rozmiar zbioru to: 268725, 336. Poni¿ej przedstawiona jest tabela zawieraj¹ca podstawowe statystyki.


```r
dim(data_final)
```

```
## [1] 268725    336
```

```r
knitr::kable(summary(data_final, digits = 2)) %>% 
  kable_styling(full_width = F) %>% 
  scroll_box(width = "100%")
```

<div style="border: 1px solid #ddd; padding: 5px; overflow-x: scroll; width:100%; "><table class="table" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:left;">   res_name </th>
   <th style="text-align:left;"> local_res_atom_non_h_count </th>
   <th style="text-align:left;"> local_res_atom_non_h_electron_sum </th>
   <th style="text-align:left;"> dict_atom_non_h_count </th>
   <th style="text-align:left;"> dict_atom_non_h_electron_sum </th>
   <th style="text-align:left;">  local_volume </th>
   <th style="text-align:left;"> local_electrons </th>
   <th style="text-align:left;">   local_mean </th>
   <th style="text-align:left;">   local_std </th>
   <th style="text-align:left;">   local_min </th>
   <th style="text-align:left;">   local_max </th>
   <th style="text-align:left;"> local_skewness </th>
   <th style="text-align:left;"> part_00_shape_segments_count </th>
   <th style="text-align:left;"> part_00_density_segments_count </th>
   <th style="text-align:left;"> part_00_volume </th>
   <th style="text-align:left;"> part_00_electrons </th>
   <th style="text-align:left;">  part_00_mean </th>
   <th style="text-align:left;">  part_00_std </th>
   <th style="text-align:left;">  part_00_max </th>
   <th style="text-align:left;"> part_00_max_over_std </th>
   <th style="text-align:left;"> part_00_skewness </th>
   <th style="text-align:left;"> part_00_parts </th>
   <th style="text-align:left;"> part_00_shape_O3 </th>
   <th style="text-align:left;"> part_00_shape_O4 </th>
   <th style="text-align:left;"> part_00_shape_O5 </th>
   <th style="text-align:left;"> part_00_shape_FL </th>
   <th style="text-align:left;"> part_00_shape_O3_norm </th>
   <th style="text-align:left;"> part_00_shape_O4_norm </th>
   <th style="text-align:left;"> part_00_shape_O5_norm </th>
   <th style="text-align:left;"> part_00_shape_FL_norm </th>
   <th style="text-align:left;"> part_00_shape_I1 </th>
   <th style="text-align:left;"> part_00_shape_I2 </th>
   <th style="text-align:left;"> part_00_shape_I3 </th>
   <th style="text-align:left;"> part_00_shape_I4 </th>
   <th style="text-align:left;"> part_00_shape_I5 </th>
   <th style="text-align:left;"> part_00_shape_I6 </th>
   <th style="text-align:left;"> part_00_shape_I1_norm </th>
   <th style="text-align:left;"> part_00_shape_I2_norm </th>
   <th style="text-align:left;"> part_00_shape_I3_norm </th>
   <th style="text-align:left;"> part_00_shape_I4_norm </th>
   <th style="text-align:left;"> part_00_shape_I5_norm </th>
   <th style="text-align:left;"> part_00_shape_I6_norm </th>
   <th style="text-align:left;"> part_00_shape_M000 </th>
   <th style="text-align:left;"> part_00_shape_CI </th>
   <th style="text-align:left;"> part_00_shape_E3_E1 </th>
   <th style="text-align:left;"> part_00_shape_E2_E1 </th>
   <th style="text-align:left;"> part_00_shape_E3_E2 </th>
   <th style="text-align:left;"> part_00_shape_sqrt_E1 </th>
   <th style="text-align:left;"> part_00_shape_sqrt_E2 </th>
   <th style="text-align:left;"> part_00_shape_sqrt_E3 </th>
   <th style="text-align:left;"> part_00_density_O3 </th>
   <th style="text-align:left;"> part_00_density_O4 </th>
   <th style="text-align:left;"> part_00_density_O5 </th>
   <th style="text-align:left;"> part_00_density_FL </th>
   <th style="text-align:left;"> part_00_density_O3_norm </th>
   <th style="text-align:left;"> part_00_density_O4_norm </th>
   <th style="text-align:left;"> part_00_density_O5_norm </th>
   <th style="text-align:left;"> part_00_density_FL_norm </th>
   <th style="text-align:left;"> part_00_density_I1 </th>
   <th style="text-align:left;"> part_00_density_I2 </th>
   <th style="text-align:left;"> part_00_density_I3 </th>
   <th style="text-align:left;"> part_00_density_I4 </th>
   <th style="text-align:left;"> part_00_density_I5 </th>
   <th style="text-align:left;"> part_00_density_I6 </th>
   <th style="text-align:left;"> part_00_density_I1_norm </th>
   <th style="text-align:left;"> part_00_density_I2_norm </th>
   <th style="text-align:left;"> part_00_density_I3_norm </th>
   <th style="text-align:left;"> part_00_density_I4_norm </th>
   <th style="text-align:left;"> part_00_density_I5_norm </th>
   <th style="text-align:left;"> part_00_density_I6_norm </th>
   <th style="text-align:left;"> part_00_density_M000 </th>
   <th style="text-align:left;"> part_00_density_CI </th>
   <th style="text-align:left;"> part_00_density_E3_E1 </th>
   <th style="text-align:left;"> part_00_density_E2_E1 </th>
   <th style="text-align:left;"> part_00_density_E3_E2 </th>
   <th style="text-align:left;"> part_00_density_sqrt_E1 </th>
   <th style="text-align:left;"> part_00_density_sqrt_E2 </th>
   <th style="text-align:left;"> part_00_density_sqrt_E3 </th>
   <th style="text-align:left;"> part_00_shape_Z_7_3 </th>
   <th style="text-align:left;"> part_00_shape_Z_0_0 </th>
   <th style="text-align:left;"> part_00_shape_Z_7_0 </th>
   <th style="text-align:left;"> part_00_shape_Z_7_1 </th>
   <th style="text-align:left;"> part_00_shape_Z_3_0 </th>
   <th style="text-align:left;"> part_00_shape_Z_5_2 </th>
   <th style="text-align:left;"> part_00_shape_Z_6_1 </th>
   <th style="text-align:left;"> part_00_shape_Z_3_1 </th>
   <th style="text-align:left;"> part_00_shape_Z_6_0 </th>
   <th style="text-align:left;"> part_00_shape_Z_2_1 </th>
   <th style="text-align:left;"> part_00_shape_Z_6_3 </th>
   <th style="text-align:left;"> part_00_shape_Z_2_0 </th>
   <th style="text-align:left;"> part_00_shape_Z_6_2 </th>
   <th style="text-align:left;"> part_00_shape_Z_5_0 </th>
   <th style="text-align:left;"> part_00_shape_Z_5_1 </th>
   <th style="text-align:left;"> part_00_shape_Z_4_2 </th>
   <th style="text-align:left;"> part_00_shape_Z_1_0 </th>
   <th style="text-align:left;"> part_00_shape_Z_4_1 </th>
   <th style="text-align:left;"> part_00_shape_Z_7_2 </th>
   <th style="text-align:left;"> part_00_shape_Z_4_0 </th>
   <th style="text-align:left;"> part_00_density_Z_7_3 </th>
   <th style="text-align:left;"> part_00_density_Z_0_0 </th>
   <th style="text-align:left;"> part_00_density_Z_7_0 </th>
   <th style="text-align:left;"> part_00_density_Z_7_1 </th>
   <th style="text-align:left;"> part_00_density_Z_3_0 </th>
   <th style="text-align:left;"> part_00_density_Z_5_2 </th>
   <th style="text-align:left;"> part_00_density_Z_6_1 </th>
   <th style="text-align:left;"> part_00_density_Z_3_1 </th>
   <th style="text-align:left;"> part_00_density_Z_6_0 </th>
   <th style="text-align:left;"> part_00_density_Z_2_1 </th>
   <th style="text-align:left;"> part_00_density_Z_6_3 </th>
   <th style="text-align:left;"> part_00_density_Z_2_0 </th>
   <th style="text-align:left;"> part_00_density_Z_6_2 </th>
   <th style="text-align:left;"> part_00_density_Z_5_0 </th>
   <th style="text-align:left;"> part_00_density_Z_5_1 </th>
   <th style="text-align:left;"> part_00_density_Z_4_2 </th>
   <th style="text-align:left;"> part_00_density_Z_1_0 </th>
   <th style="text-align:left;"> part_00_density_Z_4_1 </th>
   <th style="text-align:left;"> part_00_density_Z_7_2 </th>
   <th style="text-align:left;"> part_00_density_Z_4_0 </th>
   <th style="text-align:left;"> part_01_shape_segments_count </th>
   <th style="text-align:left;"> part_01_density_segments_count </th>
   <th style="text-align:left;"> part_01_volume </th>
   <th style="text-align:left;"> part_01_electrons </th>
   <th style="text-align:left;">  part_01_mean </th>
   <th style="text-align:left;">  part_01_std </th>
   <th style="text-align:left;">  part_01_max </th>
   <th style="text-align:left;"> part_01_max_over_std </th>
   <th style="text-align:left;"> part_01_skewness </th>
   <th style="text-align:left;"> part_01_parts </th>
   <th style="text-align:left;"> part_01_shape_O3 </th>
   <th style="text-align:left;"> part_01_shape_O4 </th>
   <th style="text-align:left;"> part_01_shape_O5 </th>
   <th style="text-align:left;"> part_01_shape_FL </th>
   <th style="text-align:left;"> part_01_shape_O3_norm </th>
   <th style="text-align:left;"> part_01_shape_O4_norm </th>
   <th style="text-align:left;"> part_01_shape_O5_norm </th>
   <th style="text-align:left;"> part_01_shape_FL_norm </th>
   <th style="text-align:left;"> part_01_shape_I1 </th>
   <th style="text-align:left;"> part_01_shape_I2 </th>
   <th style="text-align:left;"> part_01_shape_I3 </th>
   <th style="text-align:left;"> part_01_shape_I4 </th>
   <th style="text-align:left;"> part_01_shape_I5 </th>
   <th style="text-align:left;"> part_01_shape_I6 </th>
   <th style="text-align:left;"> part_01_shape_I1_norm </th>
   <th style="text-align:left;"> part_01_shape_I2_norm </th>
   <th style="text-align:left;"> part_01_shape_I3_norm </th>
   <th style="text-align:left;"> part_01_shape_I4_norm </th>
   <th style="text-align:left;"> part_01_shape_I5_norm </th>
   <th style="text-align:left;"> part_01_shape_I6_norm </th>
   <th style="text-align:left;"> part_01_shape_M000 </th>
   <th style="text-align:left;"> part_01_shape_CI </th>
   <th style="text-align:left;"> part_01_shape_E3_E1 </th>
   <th style="text-align:left;"> part_01_shape_E2_E1 </th>
   <th style="text-align:left;"> part_01_shape_E3_E2 </th>
   <th style="text-align:left;"> part_01_shape_sqrt_E1 </th>
   <th style="text-align:left;"> part_01_shape_sqrt_E2 </th>
   <th style="text-align:left;"> part_01_shape_sqrt_E3 </th>
   <th style="text-align:left;"> part_01_density_O3 </th>
   <th style="text-align:left;"> part_01_density_O4 </th>
   <th style="text-align:left;"> part_01_density_O5 </th>
   <th style="text-align:left;"> part_01_density_FL </th>
   <th style="text-align:left;"> part_01_density_O3_norm </th>
   <th style="text-align:left;"> part_01_density_O4_norm </th>
   <th style="text-align:left;"> part_01_density_O5_norm </th>
   <th style="text-align:left;"> part_01_density_FL_norm </th>
   <th style="text-align:left;"> part_01_density_I1 </th>
   <th style="text-align:left;"> part_01_density_I2 </th>
   <th style="text-align:left;"> part_01_density_I3 </th>
   <th style="text-align:left;"> part_01_density_I4 </th>
   <th style="text-align:left;"> part_01_density_I5 </th>
   <th style="text-align:left;"> part_01_density_I6 </th>
   <th style="text-align:left;"> part_01_density_I1_norm </th>
   <th style="text-align:left;"> part_01_density_I2_norm </th>
   <th style="text-align:left;"> part_01_density_I3_norm </th>
   <th style="text-align:left;"> part_01_density_I4_norm </th>
   <th style="text-align:left;"> part_01_density_I5_norm </th>
   <th style="text-align:left;"> part_01_density_I6_norm </th>
   <th style="text-align:left;"> part_01_density_M000 </th>
   <th style="text-align:left;"> part_01_density_CI </th>
   <th style="text-align:left;"> part_01_density_E3_E1 </th>
   <th style="text-align:left;"> part_01_density_E2_E1 </th>
   <th style="text-align:left;"> part_01_density_E3_E2 </th>
   <th style="text-align:left;"> part_01_density_sqrt_E1 </th>
   <th style="text-align:left;"> part_01_density_sqrt_E2 </th>
   <th style="text-align:left;"> part_01_density_sqrt_E3 </th>
   <th style="text-align:left;"> part_01_shape_Z_7_3 </th>
   <th style="text-align:left;"> part_01_shape_Z_0_0 </th>
   <th style="text-align:left;"> part_01_shape_Z_7_0 </th>
   <th style="text-align:left;"> part_01_shape_Z_7_1 </th>
   <th style="text-align:left;"> part_01_shape_Z_3_0 </th>
   <th style="text-align:left;"> part_01_shape_Z_5_2 </th>
   <th style="text-align:left;"> part_01_shape_Z_6_1 </th>
   <th style="text-align:left;"> part_01_shape_Z_3_1 </th>
   <th style="text-align:left;"> part_01_shape_Z_6_0 </th>
   <th style="text-align:left;"> part_01_shape_Z_2_1 </th>
   <th style="text-align:left;"> part_01_shape_Z_6_3 </th>
   <th style="text-align:left;"> part_01_shape_Z_2_0 </th>
   <th style="text-align:left;"> part_01_shape_Z_6_2 </th>
   <th style="text-align:left;"> part_01_shape_Z_5_0 </th>
   <th style="text-align:left;"> part_01_shape_Z_5_1 </th>
   <th style="text-align:left;"> part_01_shape_Z_4_2 </th>
   <th style="text-align:left;"> part_01_shape_Z_1_0 </th>
   <th style="text-align:left;"> part_01_shape_Z_4_1 </th>
   <th style="text-align:left;"> part_01_shape_Z_7_2 </th>
   <th style="text-align:left;"> part_01_shape_Z_4_0 </th>
   <th style="text-align:left;"> part_01_density_Z_7_3 </th>
   <th style="text-align:left;"> part_01_density_Z_0_0 </th>
   <th style="text-align:left;"> part_01_density_Z_7_0 </th>
   <th style="text-align:left;"> part_01_density_Z_7_1 </th>
   <th style="text-align:left;"> part_01_density_Z_3_0 </th>
   <th style="text-align:left;"> part_01_density_Z_5_2 </th>
   <th style="text-align:left;"> part_01_density_Z_6_1 </th>
   <th style="text-align:left;"> part_01_density_Z_3_1 </th>
   <th style="text-align:left;"> part_01_density_Z_6_0 </th>
   <th style="text-align:left;"> part_01_density_Z_2_1 </th>
   <th style="text-align:left;"> part_01_density_Z_6_3 </th>
   <th style="text-align:left;"> part_01_density_Z_2_0 </th>
   <th style="text-align:left;"> part_01_density_Z_6_2 </th>
   <th style="text-align:left;"> part_01_density_Z_5_0 </th>
   <th style="text-align:left;"> part_01_density_Z_5_1 </th>
   <th style="text-align:left;"> part_01_density_Z_4_2 </th>
   <th style="text-align:left;"> part_01_density_Z_1_0 </th>
   <th style="text-align:left;"> part_01_density_Z_4_1 </th>
   <th style="text-align:left;"> part_01_density_Z_7_2 </th>
   <th style="text-align:left;"> part_01_density_Z_4_0 </th>
   <th style="text-align:left;"> part_02_shape_segments_count </th>
   <th style="text-align:left;"> part_02_density_segments_count </th>
   <th style="text-align:left;"> part_02_volume </th>
   <th style="text-align:left;"> part_02_electrons </th>
   <th style="text-align:left;">  part_02_mean </th>
   <th style="text-align:left;">  part_02_std </th>
   <th style="text-align:left;">  part_02_max </th>
   <th style="text-align:left;"> part_02_max_over_std </th>
   <th style="text-align:left;"> part_02_skewness </th>
   <th style="text-align:left;"> part_02_parts </th>
   <th style="text-align:left;"> part_02_shape_O3 </th>
   <th style="text-align:left;"> part_02_shape_O4 </th>
   <th style="text-align:left;"> part_02_shape_O5 </th>
   <th style="text-align:left;"> part_02_shape_FL </th>
   <th style="text-align:left;"> part_02_shape_O3_norm </th>
   <th style="text-align:left;"> part_02_shape_O4_norm </th>
   <th style="text-align:left;"> part_02_shape_O5_norm </th>
   <th style="text-align:left;"> part_02_shape_FL_norm </th>
   <th style="text-align:left;"> part_02_shape_I1 </th>
   <th style="text-align:left;"> part_02_shape_I2 </th>
   <th style="text-align:left;"> part_02_shape_I3 </th>
   <th style="text-align:left;"> part_02_shape_I4 </th>
   <th style="text-align:left;"> part_02_shape_I5 </th>
   <th style="text-align:left;"> part_02_shape_I6 </th>
   <th style="text-align:left;"> part_02_shape_I1_norm </th>
   <th style="text-align:left;"> part_02_shape_I2_norm </th>
   <th style="text-align:left;"> part_02_shape_I3_norm </th>
   <th style="text-align:left;"> part_02_shape_I4_norm </th>
   <th style="text-align:left;"> part_02_shape_I5_norm </th>
   <th style="text-align:left;"> part_02_shape_I6_norm </th>
   <th style="text-align:left;"> part_02_shape_M000 </th>
   <th style="text-align:left;"> part_02_shape_CI </th>
   <th style="text-align:left;"> part_02_shape_E3_E1 </th>
   <th style="text-align:left;"> part_02_shape_E2_E1 </th>
   <th style="text-align:left;"> part_02_shape_E3_E2 </th>
   <th style="text-align:left;"> part_02_shape_sqrt_E1 </th>
   <th style="text-align:left;"> part_02_shape_sqrt_E2 </th>
   <th style="text-align:left;"> part_02_shape_sqrt_E3 </th>
   <th style="text-align:left;"> part_02_density_O3 </th>
   <th style="text-align:left;"> part_02_density_O4 </th>
   <th style="text-align:left;"> part_02_density_O5 </th>
   <th style="text-align:left;"> part_02_density_FL </th>
   <th style="text-align:left;"> part_02_density_O3_norm </th>
   <th style="text-align:left;"> part_02_density_O4_norm </th>
   <th style="text-align:left;"> part_02_density_O5_norm </th>
   <th style="text-align:left;"> part_02_density_FL_norm </th>
   <th style="text-align:left;"> part_02_density_I1 </th>
   <th style="text-align:left;"> part_02_density_I2 </th>
   <th style="text-align:left;"> part_02_density_I3 </th>
   <th style="text-align:left;"> part_02_density_I4 </th>
   <th style="text-align:left;"> part_02_density_I5 </th>
   <th style="text-align:left;"> part_02_density_I6 </th>
   <th style="text-align:left;"> part_02_density_I1_norm </th>
   <th style="text-align:left;"> part_02_density_I2_norm </th>
   <th style="text-align:left;"> part_02_density_I3_norm </th>
   <th style="text-align:left;"> part_02_density_I4_norm </th>
   <th style="text-align:left;"> part_02_density_I5_norm </th>
   <th style="text-align:left;"> part_02_density_I6_norm </th>
   <th style="text-align:left;"> part_02_density_M000 </th>
   <th style="text-align:left;"> part_02_density_CI </th>
   <th style="text-align:left;"> part_02_density_E3_E1 </th>
   <th style="text-align:left;"> part_02_density_E2_E1 </th>
   <th style="text-align:left;"> part_02_density_E3_E2 </th>
   <th style="text-align:left;"> part_02_density_sqrt_E1 </th>
   <th style="text-align:left;"> part_02_density_sqrt_E2 </th>
   <th style="text-align:left;"> part_02_density_sqrt_E3 </th>
   <th style="text-align:left;"> part_02_shape_Z_7_3 </th>
   <th style="text-align:left;"> part_02_shape_Z_0_0 </th>
   <th style="text-align:left;"> part_02_shape_Z_7_0 </th>
   <th style="text-align:left;"> part_02_shape_Z_7_1 </th>
   <th style="text-align:left;"> part_02_shape_Z_3_0 </th>
   <th style="text-align:left;"> part_02_shape_Z_5_2 </th>
   <th style="text-align:left;"> part_02_shape_Z_6_1 </th>
   <th style="text-align:left;"> part_02_shape_Z_3_1 </th>
   <th style="text-align:left;"> part_02_shape_Z_6_0 </th>
   <th style="text-align:left;"> part_02_shape_Z_2_1 </th>
   <th style="text-align:left;"> part_02_shape_Z_6_3 </th>
   <th style="text-align:left;"> part_02_shape_Z_2_0 </th>
   <th style="text-align:left;"> part_02_shape_Z_6_2 </th>
   <th style="text-align:left;"> part_02_shape_Z_5_0 </th>
   <th style="text-align:left;"> part_02_shape_Z_5_1 </th>
   <th style="text-align:left;"> part_02_shape_Z_4_2 </th>
   <th style="text-align:left;"> part_02_shape_Z_1_0 </th>
   <th style="text-align:left;"> part_02_shape_Z_4_1 </th>
   <th style="text-align:left;"> part_02_shape_Z_7_2 </th>
   <th style="text-align:left;"> part_02_shape_Z_4_0 </th>
   <th style="text-align:left;"> part_02_density_Z_7_3 </th>
   <th style="text-align:left;"> part_02_density_Z_0_0 </th>
   <th style="text-align:left;"> part_02_density_Z_7_0 </th>
   <th style="text-align:left;"> part_02_density_Z_7_1 </th>
   <th style="text-align:left;"> part_02_density_Z_3_0 </th>
   <th style="text-align:left;"> part_02_density_Z_5_2 </th>
   <th style="text-align:left;"> part_02_density_Z_6_1 </th>
   <th style="text-align:left;"> part_02_density_Z_3_1 </th>
   <th style="text-align:left;"> part_02_density_Z_6_0 </th>
   <th style="text-align:left;"> part_02_density_Z_2_1 </th>
   <th style="text-align:left;"> part_02_density_Z_6_3 </th>
   <th style="text-align:left;"> part_02_density_Z_2_0 </th>
   <th style="text-align:left;"> part_02_density_Z_6_2 </th>
   <th style="text-align:left;"> part_02_density_Z_5_0 </th>
   <th style="text-align:left;"> part_02_density_Z_5_1 </th>
   <th style="text-align:left;"> part_02_density_Z_4_2 </th>
   <th style="text-align:left;"> part_02_density_Z_1_0 </th>
   <th style="text-align:left;"> part_02_density_Z_4_1 </th>
   <th style="text-align:left;"> part_02_density_Z_7_2 </th>
   <th style="text-align:left;"> part_02_density_Z_4_0 </th>
   <th style="text-align:left;">   resolution </th>
   <th style="text-align:left;">   FoFc_mean </th>
   <th style="text-align:left;">    FoFc_std </th>
   <th style="text-align:left;"> FoFc_square_std </th>
   <th style="text-align:left;">    FoFc_min </th>
   <th style="text-align:left;">    FoFc_max </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> Length:268725 </td>
   <td style="text-align:left;"> Min.   :  1 </td>
   <td style="text-align:left;"> Min.   :   3 </td>
   <td style="text-align:left;"> Min.   :  1 </td>
   <td style="text-align:left;"> Min.   :   3 </td>
   <td style="text-align:left;"> Min.   :   93 </td>
   <td style="text-align:left;"> Min.   :  0.18 </td>
   <td style="text-align:left;"> Min.   :0.00063 </td>
   <td style="text-align:left;"> Min.   :0.0047 </td>
   <td style="text-align:left;"> Min.   :0 </td>
   <td style="text-align:left;"> Min.   : 0.04 </td>
   <td style="text-align:left;"> Min.   :0.0092 </td>
   <td style="text-align:left;"> Min.   :     1 </td>
   <td style="text-align:left;"> Min.   :     1 </td>
   <td style="text-align:left;"> Min.   :   0.82 </td>
   <td style="text-align:left;"> Min.   :  0.18 </td>
   <td style="text-align:left;"> Min.   :0.032 </td>
   <td style="text-align:left;"> Min.   :0.003 </td>
   <td style="text-align:left;"> Min.   : 0.04 </td>
   <td style="text-align:left;"> Min.   :  3.8 </td>
   <td style="text-align:left;"> Min.   :0.0023 </td>
   <td style="text-align:left;"> Min.   : 1.0 </td>
   <td style="text-align:left;"> Min.   :7.1e+02 </td>
   <td style="text-align:left;"> Min.   :1.4e+05 </td>
   <td style="text-align:left;"> Min.   :8.1e+06 </td>
   <td style="text-align:left;"> Min.   :2.0e+04 </td>
   <td style="text-align:left;"> Min.   : 0.23 </td>
   <td style="text-align:left;"> Min.   :0.018 </td>
   <td style="text-align:left;"> Min.   :0.00045 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :7.0e+03 </td>
   <td style="text-align:left;"> Min.   :9.0e+06 </td>
   <td style="text-align:left;"> Min.   :1.7e+07 </td>
   <td style="text-align:left;"> Min.   :8.4e+03 </td>
   <td style="text-align:left;"> Min.   :6.0e+00 </td>
   <td style="text-align:left;"> Min.   :2.2e+06 </td>
   <td style="text-align:left;"> Min.   :6.4e-02 </td>
   <td style="text-align:left;"> Min.   :1.1e-03 </td>
   <td style="text-align:left;"> Min.   :      0 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :   102 </td>
   <td style="text-align:left;"> Min.   :-1.3e+02 </td>
   <td style="text-align:left;"> Min.   :0.00009 </td>
   <td style="text-align:left;"> Min.   :0.00018 </td>
   <td style="text-align:left;"> Min.   :0.011 </td>
   <td style="text-align:left;"> Min.   :  1.9 </td>
   <td style="text-align:left;"> Min.   : 1.3 </td>
   <td style="text-align:left;"> Min.   : 0.85 </td>
   <td style="text-align:left;"> Min.   :3.1e+02 </td>
   <td style="text-align:left;"> Min.   :3.1e+04 </td>
   <td style="text-align:left;"> Min.   :9.2e+05 </td>
   <td style="text-align:left;"> Min.   :3.5e+03 </td>
   <td style="text-align:left;"> Min.   :  0.036 </td>
   <td style="text-align:left;"> Min.   :0.00042 </td>
   <td style="text-align:left;"> Min.   :1.7e-06 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :4.2e+03 </td>
   <td style="text-align:left;"> Min.   :4.6e+06 </td>
   <td style="text-align:left;"> Min.   :3.8e+06 </td>
   <td style="text-align:left;"> Min.   :1.4e+03 </td>
   <td style="text-align:left;"> Min.   :2.0e+00 </td>
   <td style="text-align:left;"> Min.   :4.4e+05 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :   23 </td>
   <td style="text-align:left;"> Min.   :-1.6e+02 </td>
   <td style="text-align:left;"> Min.   :8.1e-05 </td>
   <td style="text-align:left;"> Min.   :0.00016 </td>
   <td style="text-align:left;"> Min.   :0.013 </td>
   <td style="text-align:left;"> Min.   :  1.5 </td>
   <td style="text-align:left;"> Min.   : 1.2 </td>
   <td style="text-align:left;"> Min.   : 0.83 </td>
   <td style="text-align:left;"> Min.   :  7.4 </td>
   <td style="text-align:left;"> Min.   :  4.9 </td>
   <td style="text-align:left;"> Min.   :  0.85 </td>
   <td style="text-align:left;"> Min.   :  3.7 </td>
   <td style="text-align:left;"> Min.   :  0.85 </td>
   <td style="text-align:left;"> Min.   :  5.5 </td>
   <td style="text-align:left;"> Min.   :  2.2 </td>
   <td style="text-align:left;"> Min.   :  3.2 </td>
   <td style="text-align:left;"> Min.   :  0.024 </td>
   <td style="text-align:left;"> Min.   :  2.7 </td>
   <td style="text-align:left;"> Min.   :  5.5 </td>
   <td style="text-align:left;"> Min.   :  1.3 </td>
   <td style="text-align:left;"> Min.   :  4 </td>
   <td style="text-align:left;"> Min.   :  0.88 </td>
   <td style="text-align:left;"> Min.   :  3.9 </td>
   <td style="text-align:left;"> Min.   :  5.2 </td>
   <td style="text-align:left;"> Min.   :0.74 </td>
   <td style="text-align:left;"> Min.   :  3.3 </td>
   <td style="text-align:left;"> Min.   :  6.4 </td>
   <td style="text-align:left;"> Min.   :  0.03 </td>
   <td style="text-align:left;"> Min.   :  5.8 </td>
   <td style="text-align:left;"> Min.   :  2.3 </td>
   <td style="text-align:left;"> Min.   :  0.98 </td>
   <td style="text-align:left;"> Min.   :  3.1 </td>
   <td style="text-align:left;"> Min.   : 0.71 </td>
   <td style="text-align:left;"> Min.   :  3.9 </td>
   <td style="text-align:left;"> Min.   :  1.3 </td>
   <td style="text-align:left;"> Min.   :  2.0 </td>
   <td style="text-align:left;"> Min.   :7.1e-03 </td>
   <td style="text-align:left;"> Min.   :  2.9 </td>
   <td style="text-align:left;"> Min.   :  2.2 </td>
   <td style="text-align:left;"> Min.   :  2.4 </td>
   <td style="text-align:left;"> Min.   :  1.8 </td>
   <td style="text-align:left;"> Min.   :  0.87 </td>
   <td style="text-align:left;"> Min.   :  3.1 </td>
   <td style="text-align:left;"> Min.   :  2 </td>
   <td style="text-align:left;"> Min.   :0.68 </td>
   <td style="text-align:left;"> Min.   :  0.92 </td>
   <td style="text-align:left;"> Min.   :  5.1 </td>
   <td style="text-align:left;"> Min.   :7.4e-03 </td>
   <td style="text-align:left;"> Min.   :    1 </td>
   <td style="text-align:left;"> Min.   :    1 </td>
   <td style="text-align:left;"> Min.   :   0.54 </td>
   <td style="text-align:left;"> Min.   :  0.096 </td>
   <td style="text-align:left;"> Min.   :0.036 </td>
   <td style="text-align:left;"> Min.   :0.0017 </td>
   <td style="text-align:left;"> Min.   : 0.04 </td>
   <td style="text-align:left;"> Min.   :  3.8 </td>
   <td style="text-align:left;"> Min.   :0.0011 </td>
   <td style="text-align:left;"> Min.   : 1.0 </td>
   <td style="text-align:left;"> Min.   :3.0e+02 </td>
   <td style="text-align:left;"> Min.   :2.6e+04 </td>
   <td style="text-align:left;"> Min.   :6.9e+05 </td>
   <td style="text-align:left;"> Min.   :2.4e+03 </td>
   <td style="text-align:left;"> Min.   : 0.23 </td>
   <td style="text-align:left;"> Min.   : 0.018 </td>
   <td style="text-align:left;"> Min.   :0.00045 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :1.7e+03 </td>
   <td style="text-align:left;"> Min.   :5.9e+05 </td>
   <td style="text-align:left;"> Min.   :1.0e+06 </td>
   <td style="text-align:left;"> Min.   :1.0e+03 </td>
   <td style="text-align:left;"> Min.   :3.0e+00 </td>
   <td style="text-align:left;"> Min.   :2.2e+05 </td>
   <td style="text-align:left;"> Min.   :6.3e-02 </td>
   <td style="text-align:left;"> Min.   :1.1e-03 </td>
   <td style="text-align:left;"> Min.   :      0 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :  0.000 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :    67 </td>
   <td style="text-align:left;"> Min.   :-1.4e+02 </td>
   <td style="text-align:left;"> Min.   :7.3e-05 </td>
   <td style="text-align:left;"> Min.   :0.00015 </td>
   <td style="text-align:left;"> Min.   :0.012 </td>
   <td style="text-align:left;"> Min.   :  1.4 </td>
   <td style="text-align:left;"> Min.   : 0.96 </td>
   <td style="text-align:left;"> Min.   : 0.69 </td>
   <td style="text-align:left;"> Min.   :1.1e+02 </td>
   <td style="text-align:left;"> Min.   :3.9e+03 </td>
   <td style="text-align:left;"> Min.   :4.1e+04 </td>
   <td style="text-align:left;"> Min.   :4.4e+02 </td>
   <td style="text-align:left;"> Min.   :  0.035 </td>
   <td style="text-align:left;"> Min.   :4.2e-04 </td>
   <td style="text-align:left;"> Min.   :1.6e-06 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :9.2e+02 </td>
   <td style="text-align:left;"> Min.   :2.2e+05 </td>
   <td style="text-align:left;"> Min.   :1.8e+05 </td>
   <td style="text-align:left;"> Min.   :2.1e+02 </td>
   <td style="text-align:left;"> Min.   :1.0e+00 </td>
   <td style="text-align:left;"> Min.   :3.4e+04 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :   12 </td>
   <td style="text-align:left;"> Min.   :-1.6e+02 </td>
   <td style="text-align:left;"> Min.   :6.7e-05 </td>
   <td style="text-align:left;"> Min.   :0.00013 </td>
   <td style="text-align:left;"> Min.   :0.012 </td>
   <td style="text-align:left;"> Min.   :  1.4 </td>
   <td style="text-align:left;"> Min.   : 0.94 </td>
   <td style="text-align:left;"> Min.   : 0.69 </td>
   <td style="text-align:left;"> Min.   :  6.2 </td>
   <td style="text-align:left;"> Min.   :  4 </td>
   <td style="text-align:left;"> Min.   :  0.71 </td>
   <td style="text-align:left;"> Min.   :  3.9 </td>
   <td style="text-align:left;"> Min.   :  0.63 </td>
   <td style="text-align:left;"> Min.   :  4.6 </td>
   <td style="text-align:left;"> Min.   :  1.3 </td>
   <td style="text-align:left;"> Min.   :  3.1 </td>
   <td style="text-align:left;"> Min.   :  0.018 </td>
   <td style="text-align:left;"> Min.   :  3.1 </td>
   <td style="text-align:left;"> Min.   :  4.2 </td>
   <td style="text-align:left;"> Min.   :  0.37 </td>
   <td style="text-align:left;"> Min.   :  3.1 </td>
   <td style="text-align:left;"> Min.   :  0.94 </td>
   <td style="text-align:left;"> Min.   :  3.2 </td>
   <td style="text-align:left;"> Min.   :  3.6 </td>
   <td style="text-align:left;"> Min.   :0.7 </td>
   <td style="text-align:left;"> Min.   :  2 </td>
   <td style="text-align:left;"> Min.   :  5 </td>
   <td style="text-align:left;"> Min.   :  0.021 </td>
   <td style="text-align:left;"> Min.   :  5.0 </td>
   <td style="text-align:left;"> Min.   :  1.7 </td>
   <td style="text-align:left;"> Min.   :  1.3 </td>
   <td style="text-align:left;"> Min.   :  3.1 </td>
   <td style="text-align:left;"> Min.   : 0.44 </td>
   <td style="text-align:left;"> Min.   :  3.7 </td>
   <td style="text-align:left;"> Min.   :  0.74 </td>
   <td style="text-align:left;"> Min.   :  2.4 </td>
   <td style="text-align:left;"> Min.   :  0.012 </td>
   <td style="text-align:left;"> Min.   :  1.6 </td>
   <td style="text-align:left;"> Min.   :  1.5 </td>
   <td style="text-align:left;"> Min.   :  0.75 </td>
   <td style="text-align:left;"> Min.   :  1.1 </td>
   <td style="text-align:left;"> Min.   :  0.9 </td>
   <td style="text-align:left;"> Min.   :  2.3 </td>
   <td style="text-align:left;"> Min.   :  1.6 </td>
   <td style="text-align:left;"> Min.   :0.62 </td>
   <td style="text-align:left;"> Min.   :  1 </td>
   <td style="text-align:left;"> Min.   :  4.4 </td>
   <td style="text-align:left;"> Min.   :5.1e-03 </td>
   <td style="text-align:left;"> Min.   :    0 </td>
   <td style="text-align:left;"> Min.   :    0 </td>
   <td style="text-align:left;"> Min.   :   0.26 </td>
   <td style="text-align:left;"> Min.   :  0.028 </td>
   <td style="text-align:left;"> Min.   :0.039 </td>
   <td style="text-align:left;"> Min.   :0.00025 </td>
   <td style="text-align:left;"> Min.   : 0.04 </td>
   <td style="text-align:left;"> Min.   :  3.8 </td>
   <td style="text-align:left;"> Min.   :0.00017 </td>
   <td style="text-align:left;"> Min.   : 1.0 </td>
   <td style="text-align:left;"> Min.   :7.2e+01 </td>
   <td style="text-align:left;"> Min.   :1.7e+03 </td>
   <td style="text-align:left;"> Min.   :1.2e+04 </td>
   <td style="text-align:left;"> Min.   :-6.1e+01 </td>
   <td style="text-align:left;"> Min.   : 0.22 </td>
   <td style="text-align:left;"> Min.   : 0.017 </td>
   <td style="text-align:left;"> Min.   :0.00037 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :1.9e+02 </td>
   <td style="text-align:left;"> Min.   :9.3e+03 </td>
   <td style="text-align:left;"> Min.   :7.1e+03 </td>
   <td style="text-align:left;"> Min.   :-2.2e+01 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :4.5e+03 </td>
   <td style="text-align:left;"> Min.   :5.7e-02 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :    32 </td>
   <td style="text-align:left;"> Min.   :-1.5e+02 </td>
   <td style="text-align:left;"> Min.   :3.1e-05 </td>
   <td style="text-align:left;"> Min.   :4.9e-05 </td>
   <td style="text-align:left;"> Min.   :0.01 </td>
   <td style="text-align:left;"> Min.   :  0.87 </td>
   <td style="text-align:left;"> Min.   : 0.58 </td>
   <td style="text-align:left;"> Min.   : 0.42 </td>
   <td style="text-align:left;"> Min.   :1.0e+01 </td>
   <td style="text-align:left;"> Min.   :2.6e+01 </td>
   <td style="text-align:left;"> Min.   :2.0e+01 </td>
   <td style="text-align:left;"> Min.   :-1.5e+01 </td>
   <td style="text-align:left;"> Min.   :  0.035 </td>
   <td style="text-align:left;"> Min.   :4.2e-04 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :2.7e+01 </td>
   <td style="text-align:left;"> Min.   :1.8e+02 </td>
   <td style="text-align:left;"> Min.   :1.6e+02 </td>
   <td style="text-align:left;"> Min.   :-6.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :8.8e+01 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :0.0e+00 </td>
   <td style="text-align:left;"> Min.   :    3.5 </td>
   <td style="text-align:left;"> Min.   :-1.7e+02 </td>
   <td style="text-align:left;"> Min.   :3.1e-05 </td>
   <td style="text-align:left;"> Min.   :4.9e-05 </td>
   <td style="text-align:left;"> Min.   :0.0099 </td>
   <td style="text-align:left;"> Min.   :  0.87 </td>
   <td style="text-align:left;"> Min.   : 0.58 </td>
   <td style="text-align:left;"> Min.   : 0.42 </td>
   <td style="text-align:left;"> Min.   :  5.8 </td>
   <td style="text-align:left;"> Min.   :  2.8 </td>
   <td style="text-align:left;"> Min.   :  0.91 </td>
   <td style="text-align:left;"> Min.   :  3.8 </td>
   <td style="text-align:left;"> Min.   :  0.66 </td>
   <td style="text-align:left;"> Min.   :  4 </td>
   <td style="text-align:left;"> Min.   :  0.97 </td>
   <td style="text-align:left;"> Min.   :  2.6 </td>
   <td style="text-align:left;"> Min.   :2.4e-03 </td>
   <td style="text-align:left;"> Min.   :  1.6 </td>
   <td style="text-align:left;"> Min.   :  3.2 </td>
   <td style="text-align:left;"> Min.   :  0.061 </td>
   <td style="text-align:left;"> Min.   :  2.2 </td>
   <td style="text-align:left;"> Min.   :  0.88 </td>
   <td style="text-align:left;"> Min.   :  2.7 </td>
   <td style="text-align:left;"> Min.   :  2.2 </td>
   <td style="text-align:left;"> Min.   :0.67 </td>
   <td style="text-align:left;"> Min.   :  0.88 </td>
   <td style="text-align:left;"> Min.   :  4.9 </td>
   <td style="text-align:left;"> Min.   :9.5e-03 </td>
   <td style="text-align:left;"> Min.   :  3.2 </td>
   <td style="text-align:left;"> Min.   :  0.92 </td>
   <td style="text-align:left;"> Min.   :  1.2 </td>
   <td style="text-align:left;"> Min.   :  2.0 </td>
   <td style="text-align:left;"> Min.   : 0.53 </td>
   <td style="text-align:left;"> Min.   :  2.3 </td>
   <td style="text-align:left;"> Min.   :  0.46 </td>
   <td style="text-align:left;"> Min.   :  2.0 </td>
   <td style="text-align:left;"> Min.   :5.7e-03 </td>
   <td style="text-align:left;"> Min.   :  0.78 </td>
   <td style="text-align:left;"> Min.   :  1.1 </td>
   <td style="text-align:left;"> Min.   :  0.029 </td>
   <td style="text-align:left;"> Min.   :  0.82 </td>
   <td style="text-align:left;"> Min.   :  0.64 </td>
   <td style="text-align:left;"> Min.   :  1.7 </td>
   <td style="text-align:left;"> Min.   :  0.78 </td>
   <td style="text-align:left;"> Min.   :0.61 </td>
   <td style="text-align:left;"> Min.   :  0.4 </td>
   <td style="text-align:left;"> Min.   :  2.6 </td>
   <td style="text-align:left;"> Min.   :6.4e-03 </td>
   <td style="text-align:left;"> Min.   :0.48 </td>
   <td style="text-align:left;"> Min.   :-1.1e-07 </td>
   <td style="text-align:left;"> Min.   :0.010 </td>
   <td style="text-align:left;"> Min.   :0.0001 </td>
   <td style="text-align:left;"> Min.   :-7.483 </td>
   <td style="text-align:left;"> Min.   : 0.048 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> Class :character </td>
   <td style="text-align:left;"> 1st Qu.:  4 </td>
   <td style="text-align:left;"> 1st Qu.:  30 </td>
   <td style="text-align:left;"> 1st Qu.:  4 </td>
   <td style="text-align:left;"> 1st Qu.:  30 </td>
   <td style="text-align:left;"> 1st Qu.:  225 </td>
   <td style="text-align:left;"> 1st Qu.:  4.07 </td>
   <td style="text-align:left;"> 1st Qu.:0.01322 </td>
   <td style="text-align:left;"> 1st Qu.:0.0732 </td>
   <td style="text-align:left;"> 1st Qu.:0 </td>
   <td style="text-align:left;"> 1st Qu.: 0.61 </td>
   <td style="text-align:left;"> 1st Qu.:0.1289 </td>
   <td style="text-align:left;"> 1st Qu.:     5 </td>
   <td style="text-align:left;"> 1st Qu.:     5 </td>
   <td style="text-align:left;"> 1st Qu.:   7.81 </td>
   <td style="text-align:left;"> 1st Qu.:  4.05 </td>
   <td style="text-align:left;"> 1st Qu.:0.381 </td>
   <td style="text-align:left;"> 1st Qu.:0.074 </td>
   <td style="text-align:left;"> 1st Qu.: 0.61 </td>
   <td style="text-align:left;"> 1st Qu.:  5.6 </td>
   <td style="text-align:left;"> 1st Qu.:0.0650 </td>
   <td style="text-align:left;"> 1st Qu.: 1.0 </td>
   <td style="text-align:left;"> 1st Qu.:3.3e+04 </td>
   <td style="text-align:left;"> 1st Qu.:2.8e+08 </td>
   <td style="text-align:left;"> 1st Qu.:6.6e+11 </td>
   <td style="text-align:left;"> 1st Qu.:1.5e+09 </td>
   <td style="text-align:left;"> 1st Qu.: 0.27 </td>
   <td style="text-align:left;"> 1st Qu.:0.022 </td>
   <td style="text-align:left;"> 1st Qu.:0.00052 </td>
   <td style="text-align:left;"> 1st Qu.:6.4e-04 </td>
   <td style="text-align:left;"> 1st Qu.:1.6e+06 </td>
   <td style="text-align:left;"> 1st Qu.:3.9e+11 </td>
   <td style="text-align:left;"> 1st Qu.:8.7e+11 </td>
   <td style="text-align:left;"> 1st Qu.:7.8e+08 </td>
   <td style="text-align:left;"> 1st Qu.:1.1e+08 </td>
   <td style="text-align:left;"> 1st Qu.:2.3e+10 </td>
   <td style="text-align:left;"> 1st Qu.:9.8e-02 </td>
   <td style="text-align:left;"> 1st Qu.:1.9e-03 </td>
   <td style="text-align:left;"> 1st Qu.:      0 </td>
   <td style="text-align:left;"> 1st Qu.:3.0e-04 </td>
   <td style="text-align:left;"> 1st Qu.:4.0e-05 </td>
   <td style="text-align:left;"> 1st Qu.:1.0e-02 </td>
   <td style="text-align:left;"> 1st Qu.:   976 </td>
   <td style="text-align:left;"> 1st Qu.:-7.8e-01 </td>
   <td style="text-align:left;"> 1st Qu.:0.08866 </td>
   <td style="text-align:left;"> 1st Qu.:0.21994 </td>
   <td style="text-align:left;"> 1st Qu.:0.372 </td>
   <td style="text-align:left;"> 1st Qu.:  4.1 </td>
   <td style="text-align:left;"> 1st Qu.: 2.7 </td>
   <td style="text-align:left;"> 1st Qu.: 2.05 </td>
   <td style="text-align:left;"> 1st Qu.:1.6e+04 </td>
   <td style="text-align:left;"> 1st Qu.:6.7e+07 </td>
   <td style="text-align:left;"> 1st Qu.:7.6e+10 </td>
   <td style="text-align:left;"> 1st Qu.:3.1e+08 </td>
   <td style="text-align:left;"> 1st Qu.:  0.368 </td>
   <td style="text-align:left;"> 1st Qu.:0.03857 </td>
   <td style="text-align:left;"> 1st Qu.:1.2e-03 </td>
   <td style="text-align:left;"> 1st Qu.:1.8e-03 </td>
   <td style="text-align:left;"> 1st Qu.:7.4e+05 </td>
   <td style="text-align:left;"> 1st Qu.:8.6e+10 </td>
   <td style="text-align:left;"> 1st Qu.:2.0e+11 </td>
   <td style="text-align:left;"> 1st Qu.:1.7e+08 </td>
   <td style="text-align:left;"> 1st Qu.:3.4e+07 </td>
   <td style="text-align:left;"> 1st Qu.:5.4e+09 </td>
   <td style="text-align:left;"> 1st Qu.:2.0e-01 </td>
   <td style="text-align:left;"> 1st Qu.:7.5e-03 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:9.0e-04 </td>
   <td style="text-align:left;"> 1st Qu.:2.0e-04 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:  506 </td>
   <td style="text-align:left;"> 1st Qu.:-8.5e-01 </td>
   <td style="text-align:left;"> 1st Qu.:8.5e-02 </td>
   <td style="text-align:left;"> 1st Qu.:0.21494 </td>
   <td style="text-align:left;"> 1st Qu.:0.371 </td>
   <td style="text-align:left;"> 1st Qu.:  3.8 </td>
   <td style="text-align:left;"> 1st Qu.: 2.6 </td>
   <td style="text-align:left;"> 1st Qu.: 1.96 </td>
   <td style="text-align:left;"> 1st Qu.: 16.6 </td>
   <td style="text-align:left;"> 1st Qu.: 15.3 </td>
   <td style="text-align:left;"> 1st Qu.:  6.59 </td>
   <td style="text-align:left;"> 1st Qu.: 10.5 </td>
   <td style="text-align:left;"> 1st Qu.:  6.34 </td>
   <td style="text-align:left;"> 1st Qu.: 15.8 </td>
   <td style="text-align:left;"> 1st Qu.: 12.5 </td>
   <td style="text-align:left;"> 1st Qu.: 12.2 </td>
   <td style="text-align:left;"> 1st Qu.:  5.774 </td>
   <td style="text-align:left;"> 1st Qu.: 20.7 </td>
   <td style="text-align:left;"> 1st Qu.: 19.6 </td>
   <td style="text-align:left;"> 1st Qu.: 15.4 </td>
   <td style="text-align:left;"> 1st Qu.: 17 </td>
   <td style="text-align:left;"> 1st Qu.:  5.79 </td>
   <td style="text-align:left;"> 1st Qu.: 12.4 </td>
   <td style="text-align:left;"> 1st Qu.: 20.3 </td>
   <td style="text-align:left;"> 1st Qu.:1.26 </td>
   <td style="text-align:left;"> 1st Qu.: 17.0 </td>
   <td style="text-align:left;"> 1st Qu.: 14.1 </td>
   <td style="text-align:left;"> 1st Qu.:  8.67 </td>
   <td style="text-align:left;"> 1st Qu.: 11.4 </td>
   <td style="text-align:left;"> 1st Qu.: 11.0 </td>
   <td style="text-align:left;"> 1st Qu.:  6.05 </td>
   <td style="text-align:left;"> 1st Qu.:  7.8 </td>
   <td style="text-align:left;"> 1st Qu.: 4.74 </td>
   <td style="text-align:left;"> 1st Qu.: 10.9 </td>
   <td style="text-align:left;"> 1st Qu.:  8.3 </td>
   <td style="text-align:left;"> 1st Qu.:  8.2 </td>
   <td style="text-align:left;"> 1st Qu.:3.9e+00 </td>
   <td style="text-align:left;"> 1st Qu.: 15.4 </td>
   <td style="text-align:left;"> 1st Qu.: 13.1 </td>
   <td style="text-align:left;"> 1st Qu.: 12.0 </td>
   <td style="text-align:left;"> 1st Qu.: 11.7 </td>
   <td style="text-align:left;"> 1st Qu.:  5.31 </td>
   <td style="text-align:left;"> 1st Qu.:  9.0 </td>
   <td style="text-align:left;"> 1st Qu.: 16 </td>
   <td style="text-align:left;"> 1st Qu.:1.25 </td>
   <td style="text-align:left;"> 1st Qu.: 13.84 </td>
   <td style="text-align:left;"> 1st Qu.: 10.0 </td>
   <td style="text-align:left;"> 1st Qu.:7.9e+00 </td>
   <td style="text-align:left;"> 1st Qu.:    3 </td>
   <td style="text-align:left;"> 1st Qu.:    3 </td>
   <td style="text-align:left;"> 1st Qu.:   5.13 </td>
   <td style="text-align:left;"> 1st Qu.:  2.841 </td>
   <td style="text-align:left;"> 1st Qu.:0.418 </td>
   <td style="text-align:left;"> 1st Qu.:0.0616 </td>
   <td style="text-align:left;"> 1st Qu.: 0.61 </td>
   <td style="text-align:left;"> 1st Qu.:  5.6 </td>
   <td style="text-align:left;"> 1st Qu.:0.0525 </td>
   <td style="text-align:left;"> 1st Qu.: 1.0 </td>
   <td style="text-align:left;"> 1st Qu.:1.7e+04 </td>
   <td style="text-align:left;"> 1st Qu.:7.1e+07 </td>
   <td style="text-align:left;"> 1st Qu.:8.2e+10 </td>
   <td style="text-align:left;"> 1st Qu.:2.1e+08 </td>
   <td style="text-align:left;"> 1st Qu.: 0.26 </td>
   <td style="text-align:left;"> 1st Qu.: 0.021 </td>
   <td style="text-align:left;"> 1st Qu.:0.00049 </td>
   <td style="text-align:left;"> 1st Qu.:4.0e-04 </td>
   <td style="text-align:left;"> 1st Qu.:6.2e+05 </td>
   <td style="text-align:left;"> 1st Qu.:5.8e+10 </td>
   <td style="text-align:left;"> 1st Qu.:1.4e+11 </td>
   <td style="text-align:left;"> 1st Qu.:1.0e+08 </td>
   <td style="text-align:left;"> 1st Qu.:1.2e+07 </td>
   <td style="text-align:left;"> 1st Qu.:4.7e+09 </td>
   <td style="text-align:left;"> 1st Qu.:8.8e-02 </td>
   <td style="text-align:left;"> 1st Qu.:1.7e-03 </td>
   <td style="text-align:left;"> 1st Qu.:      0 </td>
   <td style="text-align:left;"> 1st Qu.:2.0e-04 </td>
   <td style="text-align:left;"> 1st Qu.:  0.000 </td>
   <td style="text-align:left;"> 1st Qu.:1.0e-02 </td>
   <td style="text-align:left;"> 1st Qu.:   641 </td>
   <td style="text-align:left;"> 1st Qu.:-6.3e-01 </td>
   <td style="text-align:left;"> 1st Qu.:8.0e-02 </td>
   <td style="text-align:left;"> 1st Qu.:0.20797 </td>
   <td style="text-align:left;"> 1st Qu.:0.378 </td>
   <td style="text-align:left;"> 1st Qu.:  3.6 </td>
   <td style="text-align:left;"> 1st Qu.: 2.35 </td>
   <td style="text-align:left;"> 1st Qu.: 1.78 </td>
   <td style="text-align:left;"> 1st Qu.:9.3e+03 </td>
   <td style="text-align:left;"> 1st Qu.:2.1e+07 </td>
   <td style="text-align:left;"> 1st Qu.:1.3e+10 </td>
   <td style="text-align:left;"> 1st Qu.:5.2e+07 </td>
   <td style="text-align:left;"> 1st Qu.:  0.343 </td>
   <td style="text-align:left;"> 1st Qu.:3.4e-02 </td>
   <td style="text-align:left;"> 1st Qu.:9.8e-04 </td>
   <td style="text-align:left;"> 1st Qu.:1.0e-03 </td>
   <td style="text-align:left;"> 1st Qu.:3.3e+05 </td>
   <td style="text-align:left;"> 1st Qu.:1.6e+10 </td>
   <td style="text-align:left;"> 1st Qu.:3.8e+10 </td>
   <td style="text-align:left;"> 1st Qu.:2.7e+07 </td>
   <td style="text-align:left;"> 1st Qu.:4.5e+06 </td>
   <td style="text-align:left;"> 1st Qu.:1.4e+09 </td>
   <td style="text-align:left;"> 1st Qu.:1.7e-01 </td>
   <td style="text-align:left;"> 1st Qu.:6.0e-03 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:  355 </td>
   <td style="text-align:left;"> 1st Qu.:-6.7e-01 </td>
   <td style="text-align:left;"> 1st Qu.:7.8e-02 </td>
   <td style="text-align:left;"> 1st Qu.:0.20445 </td>
   <td style="text-align:left;"> 1st Qu.:0.379 </td>
   <td style="text-align:left;"> 1st Qu.:  3.4 </td>
   <td style="text-align:left;"> 1st Qu.: 2.24 </td>
   <td style="text-align:left;"> 1st Qu.: 1.72 </td>
   <td style="text-align:left;"> 1st Qu.: 13.1 </td>
   <td style="text-align:left;"> 1st Qu.: 12 </td>
   <td style="text-align:left;"> 1st Qu.:  6.59 </td>
   <td style="text-align:left;"> 1st Qu.:  8.7 </td>
   <td style="text-align:left;"> 1st Qu.:  4.93 </td>
   <td style="text-align:left;"> 1st Qu.: 12.0 </td>
   <td style="text-align:left;"> 1st Qu.:  9.9 </td>
   <td style="text-align:left;"> 1st Qu.:  9.6 </td>
   <td style="text-align:left;"> 1st Qu.:  4.625 </td>
   <td style="text-align:left;"> 1st Qu.: 16.5 </td>
   <td style="text-align:left;"> 1st Qu.: 15.3 </td>
   <td style="text-align:left;"> 1st Qu.: 11.98 </td>
   <td style="text-align:left;"> 1st Qu.: 13.3 </td>
   <td style="text-align:left;"> 1st Qu.:  5.41 </td>
   <td style="text-align:left;"> 1st Qu.:  9.2 </td>
   <td style="text-align:left;"> 1st Qu.: 15.6 </td>
   <td style="text-align:left;"> 1st Qu.:1.3 </td>
   <td style="text-align:left;"> 1st Qu.: 13 </td>
   <td style="text-align:left;"> 1st Qu.: 11 </td>
   <td style="text-align:left;"> 1st Qu.:  6.500 </td>
   <td style="text-align:left;"> 1st Qu.:  9.8 </td>
   <td style="text-align:left;"> 1st Qu.:  9.2 </td>
   <td style="text-align:left;"> 1st Qu.:  6.2 </td>
   <td style="text-align:left;"> 1st Qu.:  7.5 </td>
   <td style="text-align:left;"> 1st Qu.: 4.14 </td>
   <td style="text-align:left;"> 1st Qu.:  8.9 </td>
   <td style="text-align:left;"> 1st Qu.:  6.43 </td>
   <td style="text-align:left;"> 1st Qu.:  6.9 </td>
   <td style="text-align:left;"> 1st Qu.:  3.150 </td>
   <td style="text-align:left;"> 1st Qu.: 12.8 </td>
   <td style="text-align:left;"> 1st Qu.: 10.3 </td>
   <td style="text-align:left;"> 1st Qu.:  9.70 </td>
   <td style="text-align:left;"> 1st Qu.:  8.9 </td>
   <td style="text-align:left;"> 1st Qu.:  5.2 </td>
   <td style="text-align:left;"> 1st Qu.:  7.2 </td>
   <td style="text-align:left;"> 1st Qu.: 12.1 </td>
   <td style="text-align:left;"> 1st Qu.:1.28 </td>
   <td style="text-align:left;"> 1st Qu.: 10 </td>
   <td style="text-align:left;"> 1st Qu.:  8.6 </td>
   <td style="text-align:left;"> 1st Qu.:5.5e+00 </td>
   <td style="text-align:left;"> 1st Qu.:    3 </td>
   <td style="text-align:left;"> 1st Qu.:    3 </td>
   <td style="text-align:left;"> 1st Qu.:   3.10 </td>
   <td style="text-align:left;"> 1st Qu.:  1.832 </td>
   <td style="text-align:left;"> 1st Qu.:0.453 </td>
   <td style="text-align:left;"> 1st Qu.:0.04870 </td>
   <td style="text-align:left;"> 1st Qu.: 0.61 </td>
   <td style="text-align:left;"> 1st Qu.:  5.6 </td>
   <td style="text-align:left;"> 1st Qu.:0.04018 </td>
   <td style="text-align:left;"> 1st Qu.: 1.0 </td>
   <td style="text-align:left;"> 1st Qu.:7.5e+03 </td>
   <td style="text-align:left;"> 1st Qu.:1.3e+07 </td>
   <td style="text-align:left;"> 1st Qu.:6.5e+09 </td>
   <td style="text-align:left;"> 1st Qu.: 1.4e+07 </td>
   <td style="text-align:left;"> 1st Qu.: 0.25 </td>
   <td style="text-align:left;"> 1st Qu.: 0.020 </td>
   <td style="text-align:left;"> 1st Qu.:0.00047 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:1.8e+05 </td>
   <td style="text-align:left;"> 1st Qu.:5.2e+09 </td>
   <td style="text-align:left;"> 1st Qu.:1.1e+10 </td>
   <td style="text-align:left;"> 1st Qu.: 6.5e+06 </td>
   <td style="text-align:left;"> 1st Qu.:5.9e+05 </td>
   <td style="text-align:left;"> 1st Qu.:5.9e+08 </td>
   <td style="text-align:left;"> 1st Qu.:8.1e-02 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:   387 </td>
   <td style="text-align:left;"> 1st Qu.:-4.5e-01 </td>
   <td style="text-align:left;"> 1st Qu.:7.7e-02 </td>
   <td style="text-align:left;"> 1st Qu.:2.1e-01 </td>
   <td style="text-align:left;"> 1st Qu.:0.40 </td>
   <td style="text-align:left;"> 1st Qu.:  2.94 </td>
   <td style="text-align:left;"> 1st Qu.: 1.95 </td>
   <td style="text-align:left;"> 1st Qu.: 1.51 </td>
   <td style="text-align:left;"> 1st Qu.:4.4e+03 </td>
   <td style="text-align:left;"> 1st Qu.:4.4e+06 </td>
   <td style="text-align:left;"> 1st Qu.:1.3e+09 </td>
   <td style="text-align:left;"> 1st Qu.: 4.5e+06 </td>
   <td style="text-align:left;"> 1st Qu.:  0.318 </td>
   <td style="text-align:left;"> 1st Qu.:3.0e-02 </td>
   <td style="text-align:left;"> 1st Qu.:8.2e-04 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:1.1e+05 </td>
   <td style="text-align:left;"> 1st Qu.:1.7e+09 </td>
   <td style="text-align:left;"> 1st Qu.:3.7e+09 </td>
   <td style="text-align:left;"> 1st Qu.: 2.2e+06 </td>
   <td style="text-align:left;"> 1st Qu.:3.0e+05 </td>
   <td style="text-align:left;"> 1st Qu.:2.0e+08 </td>
   <td style="text-align:left;"> 1st Qu.:1.4e-01 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 1st Qu.:  229.0 </td>
   <td style="text-align:left;"> 1st Qu.:-4.7e-01 </td>
   <td style="text-align:left;"> 1st Qu.:7.6e-02 </td>
   <td style="text-align:left;"> 1st Qu.:2.0e-01 </td>
   <td style="text-align:left;"> 1st Qu.:0.4007 </td>
   <td style="text-align:left;"> 1st Qu.:  2.79 </td>
   <td style="text-align:left;"> 1st Qu.: 1.88 </td>
   <td style="text-align:left;"> 1st Qu.: 1.47 </td>
   <td style="text-align:left;"> 1st Qu.: 10.9 </td>
   <td style="text-align:left;"> 1st Qu.:  9.6 </td>
   <td style="text-align:left;"> 1st Qu.:  6.77 </td>
   <td style="text-align:left;"> 1st Qu.:  8.4 </td>
   <td style="text-align:left;"> 1st Qu.:  4.42 </td>
   <td style="text-align:left;"> 1st Qu.:  9 </td>
   <td style="text-align:left;"> 1st Qu.:  7.18 </td>
   <td style="text-align:left;"> 1st Qu.:  7.1 </td>
   <td style="text-align:left;"> 1st Qu.:3.4e+00 </td>
   <td style="text-align:left;"> 1st Qu.: 12.4 </td>
   <td style="text-align:left;"> 1st Qu.: 11.0 </td>
   <td style="text-align:left;"> 1st Qu.:  8.760 </td>
   <td style="text-align:left;"> 1st Qu.:  9.4 </td>
   <td style="text-align:left;"> 1st Qu.:  5.51 </td>
   <td style="text-align:left;"> 1st Qu.:  7.4 </td>
   <td style="text-align:left;"> 1st Qu.: 11.3 </td>
   <td style="text-align:left;"> 1st Qu.:1.34 </td>
   <td style="text-align:left;"> 1st Qu.:  8.90 </td>
   <td style="text-align:left;"> 1st Qu.:  9.7 </td>
   <td style="text-align:left;"> 1st Qu.:4.6e+00 </td>
   <td style="text-align:left;"> 1st Qu.:  9.3 </td>
   <td style="text-align:left;"> 1st Qu.:  7.39 </td>
   <td style="text-align:left;"> 1st Qu.:  6.5 </td>
   <td style="text-align:left;"> 1st Qu.:  7.7 </td>
   <td style="text-align:left;"> 1st Qu.: 4.13 </td>
   <td style="text-align:left;"> 1st Qu.:  7.7 </td>
   <td style="text-align:left;"> 1st Qu.:  5.02 </td>
   <td style="text-align:left;"> 1st Qu.:  5.8 </td>
   <td style="text-align:left;"> 1st Qu.:2.4e+00 </td>
   <td style="text-align:left;"> 1st Qu.: 10.11 </td>
   <td style="text-align:left;"> 1st Qu.:  7.7 </td>
   <td style="text-align:left;"> 1st Qu.:  7.401 </td>
   <td style="text-align:left;"> 1st Qu.:  6.60 </td>
   <td style="text-align:left;"> 1st Qu.:  5.36 </td>
   <td style="text-align:left;"> 1st Qu.:  6.7 </td>
   <td style="text-align:left;"> 1st Qu.:  8.40 </td>
   <td style="text-align:left;"> 1st Qu.:1.33 </td>
   <td style="text-align:left;"> 1st Qu.:  6.8 </td>
   <td style="text-align:left;"> 1st Qu.:  8.6 </td>
   <td style="text-align:left;"> 1st Qu.:3.3e+00 </td>
   <td style="text-align:left;"> 1st Qu.:1.80 </td>
   <td style="text-align:left;"> 1st Qu.:-4.6e-11 </td>
   <td style="text-align:left;"> 1st Qu.:0.092 </td>
   <td style="text-align:left;"> 1st Qu.:0.0085 </td>
   <td style="text-align:left;"> 1st Qu.:-0.857 </td>
   <td style="text-align:left;"> 1st Qu.: 1.186 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> Mode  :character </td>
   <td style="text-align:left;"> Median :  6 </td>
   <td style="text-align:left;"> Median :  48 </td>
   <td style="text-align:left;"> Median :  6 </td>
   <td style="text-align:left;"> Median :  48 </td>
   <td style="text-align:left;"> Median :  364 </td>
   <td style="text-align:left;"> Median :  8.61 </td>
   <td style="text-align:left;"> Median :0.01961 </td>
   <td style="text-align:left;"> Median :0.1030 </td>
   <td style="text-align:left;"> Median :0 </td>
   <td style="text-align:left;"> Median : 0.94 </td>
   <td style="text-align:left;"> Median :0.1821 </td>
   <td style="text-align:left;"> Median :    26 </td>
   <td style="text-align:left;"> Median :    26 </td>
   <td style="text-align:left;"> Median :  15.41 </td>
   <td style="text-align:left;"> Median :  8.53 </td>
   <td style="text-align:left;"> Median :0.532 </td>
   <td style="text-align:left;"> Median :0.131 </td>
   <td style="text-align:left;"> Median : 0.94 </td>
   <td style="text-align:left;"> Median :  7.6 </td>
   <td style="text-align:left;"> Median :0.1211 </td>
   <td style="text-align:left;"> Median : 1.0 </td>
   <td style="text-align:left;"> Median :1.1e+05 </td>
   <td style="text-align:left;"> Median :3.0e+09 </td>
   <td style="text-align:left;"> Median :2.3e+13 </td>
   <td style="text-align:left;"> Median :5.5e+10 </td>
   <td style="text-align:left;"> Median : 0.38 </td>
   <td style="text-align:left;"> Median :0.034 </td>
   <td style="text-align:left;"> Median :0.00079 </td>
   <td style="text-align:left;"> Median :5.7e-03 </td>
   <td style="text-align:left;"> Median :9.5e+06 </td>
   <td style="text-align:left;"> Median :1.3e+13 </td>
   <td style="text-align:left;"> Median :3.7e+13 </td>
   <td style="text-align:left;"> Median :3.1e+10 </td>
   <td style="text-align:left;"> Median :9.5e+09 </td>
   <td style="text-align:left;"> Median :4.7e+11 </td>
   <td style="text-align:left;"> Median :2.3e-01 </td>
   <td style="text-align:left;"> Median :6.8e-03 </td>
   <td style="text-align:left;"> Median :      0 </td>
   <td style="text-align:left;"> Median :3.2e-03 </td>
   <td style="text-align:left;"> Median :1.1e-03 </td>
   <td style="text-align:left;"> Median :4.0e-02 </td>
   <td style="text-align:left;"> Median :  1926 </td>
   <td style="text-align:left;"> Median : 2.1e-04 </td>
   <td style="text-align:left;"> Median :0.17421 </td>
   <td style="text-align:left;"> Median :0.39055 </td>
   <td style="text-align:left;"> Median :0.577 </td>
   <td style="text-align:left;"> Median :  6.0 </td>
   <td style="text-align:left;"> Median : 3.6 </td>
   <td style="text-align:left;"> Median : 2.67 </td>
   <td style="text-align:left;"> Median :5.4e+04 </td>
   <td style="text-align:left;"> Median :7.3e+08 </td>
   <td style="text-align:left;"> Median :2.7e+12 </td>
   <td style="text-align:left;"> Median :1.3e+10 </td>
   <td style="text-align:left;"> Median :  0.586 </td>
   <td style="text-align:left;"> Median :0.08264 </td>
   <td style="text-align:left;"> Median :2.9e-03 </td>
   <td style="text-align:left;"> Median :1.8e-02 </td>
   <td style="text-align:left;"> Median :4.2e+06 </td>
   <td style="text-align:left;"> Median :2.6e+12 </td>
   <td style="text-align:left;"> Median :7.4e+12 </td>
   <td style="text-align:left;"> Median :7.7e+09 </td>
   <td style="text-align:left;"> Median :3.1e+09 </td>
   <td style="text-align:left;"> Median :1.0e+11 </td>
   <td style="text-align:left;"> Median :5.4e-01 </td>
   <td style="text-align:left;"> Median :4.0e-02 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :1.1e-02 </td>
   <td style="text-align:left;"> Median :4.6e-03 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median : 1066 </td>
   <td style="text-align:left;"> Median : 1.1e-04 </td>
   <td style="text-align:left;"> Median :1.7e-01 </td>
   <td style="text-align:left;"> Median :0.39101 </td>
   <td style="text-align:left;"> Median :0.581 </td>
   <td style="text-align:left;"> Median :  5.7 </td>
   <td style="text-align:left;"> Median : 3.4 </td>
   <td style="text-align:left;"> Median : 2.51 </td>
   <td style="text-align:left;"> Median : 28.1 </td>
   <td style="text-align:left;"> Median : 21.4 </td>
   <td style="text-align:left;"> Median : 10.72 </td>
   <td style="text-align:left;"> Median : 18.7 </td>
   <td style="text-align:left;"> Median : 11.25 </td>
   <td style="text-align:left;"> Median : 26.0 </td>
   <td style="text-align:left;"> Median : 22.0 </td>
   <td style="text-align:left;"> Median : 19.0 </td>
   <td style="text-align:left;"> Median : 10.389 </td>
   <td style="text-align:left;"> Median : 29.9 </td>
   <td style="text-align:left;"> Median : 32.9 </td>
   <td style="text-align:left;"> Median : 22.8 </td>
   <td style="text-align:left;"> Median : 29 </td>
   <td style="text-align:left;"> Median : 12.95 </td>
   <td style="text-align:left;"> Median : 21.6 </td>
   <td style="text-align:left;"> Median : 32.9 </td>
   <td style="text-align:left;"> Median :1.40 </td>
   <td style="text-align:left;"> Median : 28.3 </td>
   <td style="text-align:left;"> Median : 24.7 </td>
   <td style="text-align:left;"> Median : 15.73 </td>
   <td style="text-align:left;"> Median : 20.8 </td>
   <td style="text-align:left;"> Median : 16.0 </td>
   <td style="text-align:left;"> Median :  9.26 </td>
   <td style="text-align:left;"> Median : 15.0 </td>
   <td style="text-align:left;"> Median : 8.27 </td>
   <td style="text-align:left;"> Median : 18.9 </td>
   <td style="text-align:left;"> Median : 18.4 </td>
   <td style="text-align:left;"> Median : 13.2 </td>
   <td style="text-align:left;"> Median :8.8e+00 </td>
   <td style="text-align:left;"> Median : 22.5 </td>
   <td style="text-align:left;"> Median : 25.1 </td>
   <td style="text-align:left;"> Median : 17.9 </td>
   <td style="text-align:left;"> Median : 23.2 </td>
   <td style="text-align:left;"> Median : 10.48 </td>
   <td style="text-align:left;"> Median : 16.3 </td>
   <td style="text-align:left;"> Median : 25 </td>
   <td style="text-align:left;"> Median :1.38 </td>
   <td style="text-align:left;"> Median : 22.15 </td>
   <td style="text-align:left;"> Median : 18.8 </td>
   <td style="text-align:left;"> Median :1.4e+01 </td>
   <td style="text-align:left;"> Median :   17 </td>
   <td style="text-align:left;"> Median :   17 </td>
   <td style="text-align:left;"> Median :  11.35 </td>
   <td style="text-align:left;"> Median :  6.868 </td>
   <td style="text-align:left;"> Median :0.582 </td>
   <td style="text-align:left;"> Median :0.1181 </td>
   <td style="text-align:left;"> Median : 0.94 </td>
   <td style="text-align:left;"> Median :  7.6 </td>
   <td style="text-align:left;"> Median :0.1070 </td>
   <td style="text-align:left;"> Median : 1.0 </td>
   <td style="text-align:left;"> Median :6.6e+04 </td>
   <td style="text-align:left;"> Median :1.1e+09 </td>
   <td style="text-align:left;"> Median :4.9e+12 </td>
   <td style="text-align:left;"> Median :1.5e+10 </td>
   <td style="text-align:left;"> Median : 0.37 </td>
   <td style="text-align:left;"> Median : 0.032 </td>
   <td style="text-align:left;"> Median :0.00073 </td>
   <td style="text-align:left;"> Median :5.0e-03 </td>
   <td style="text-align:left;"> Median :4.6e+06 </td>
   <td style="text-align:left;"> Median :3.1e+12 </td>
   <td style="text-align:left;"> Median :8.7e+12 </td>
   <td style="text-align:left;"> Median :8.2e+09 </td>
   <td style="text-align:left;"> Median :2.4e+09 </td>
   <td style="text-align:left;"> Median :1.4e+11 </td>
   <td style="text-align:left;"> Median :2.2e-01 </td>
   <td style="text-align:left;"> Median :6.0e-03 </td>
   <td style="text-align:left;"> Median :      0 </td>
   <td style="text-align:left;"> Median :2.9e-03 </td>
   <td style="text-align:left;"> Median :  0.001 </td>
   <td style="text-align:left;"> Median :4.0e-02 </td>
   <td style="text-align:left;"> Median :  1419 </td>
   <td style="text-align:left;"> Median : 8.0e-05 </td>
   <td style="text-align:left;"> Median :1.8e-01 </td>
   <td style="text-align:left;"> Median :0.39537 </td>
   <td style="text-align:left;"> Median :0.597 </td>
   <td style="text-align:left;"> Median :  5.4 </td>
   <td style="text-align:left;"> Median : 3.25 </td>
   <td style="text-align:left;"> Median : 2.42 </td>
   <td style="text-align:left;"> Median :3.7e+04 </td>
   <td style="text-align:left;"> Median :3.4e+08 </td>
   <td style="text-align:left;"> Median :8.4e+11 </td>
   <td style="text-align:left;"> Median :4.0e+09 </td>
   <td style="text-align:left;"> Median :  0.553 </td>
   <td style="text-align:left;"> Median :7.4e-02 </td>
   <td style="text-align:left;"> Median :2.6e-03 </td>
   <td style="text-align:left;"> Median :1.3e-02 </td>
   <td style="text-align:left;"> Median :2.4e+06 </td>
   <td style="text-align:left;"> Median :8.3e+11 </td>
   <td style="text-align:left;"> Median :2.2e+12 </td>
   <td style="text-align:left;"> Median :2.4e+09 </td>
   <td style="text-align:left;"> Median :9.5e+08 </td>
   <td style="text-align:left;"> Median :3.9e+10 </td>
   <td style="text-align:left;"> Median :4.8e-01 </td>
   <td style="text-align:left;"> Median :3.2e-02 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :8.0e-03 </td>
   <td style="text-align:left;"> Median :3.0e-03 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :  859 </td>
   <td style="text-align:left;"> Median : 3.0e-05 </td>
   <td style="text-align:left;"> Median :1.8e-01 </td>
   <td style="text-align:left;"> Median :0.39667 </td>
   <td style="text-align:left;"> Median :0.600 </td>
   <td style="text-align:left;"> Median :  5.1 </td>
   <td style="text-align:left;"> Median : 3.04 </td>
   <td style="text-align:left;"> Median : 2.29 </td>
   <td style="text-align:left;"> Median : 23.0 </td>
   <td style="text-align:left;"> Median : 18 </td>
   <td style="text-align:left;"> Median :  8.94 </td>
   <td style="text-align:left;"> Median : 15.1 </td>
   <td style="text-align:left;"> Median :  9.48 </td>
   <td style="text-align:left;"> Median : 21.5 </td>
   <td style="text-align:left;"> Median : 17.8 </td>
   <td style="text-align:left;"> Median : 16.1 </td>
   <td style="text-align:left;"> Median :  8.631 </td>
   <td style="text-align:left;"> Median : 25.1 </td>
   <td style="text-align:left;"> Median : 26.9 </td>
   <td style="text-align:left;"> Median : 18.90 </td>
   <td style="text-align:left;"> Median : 23.8 </td>
   <td style="text-align:left;"> Median : 10.47 </td>
   <td style="text-align:left;"> Median : 17.5 </td>
   <td style="text-align:left;"> Median : 26.9 </td>
   <td style="text-align:left;"> Median :1.5 </td>
   <td style="text-align:left;"> Median : 23 </td>
   <td style="text-align:left;"> Median : 20 </td>
   <td style="text-align:left;"> Median : 12.364 </td>
   <td style="text-align:left;"> Median : 17.6 </td>
   <td style="text-align:left;"> Median : 14.3 </td>
   <td style="text-align:left;"> Median :  7.8 </td>
   <td style="text-align:left;"> Median : 12.4 </td>
   <td style="text-align:left;"> Median : 7.29 </td>
   <td style="text-align:left;"> Median : 16.2 </td>
   <td style="text-align:left;"> Median : 14.93 </td>
   <td style="text-align:left;"> Median : 11.7 </td>
   <td style="text-align:left;"> Median :  6.802 </td>
   <td style="text-align:left;"> Median : 20.0 </td>
   <td style="text-align:left;"> Median : 21.1 </td>
   <td style="text-align:left;"> Median : 15.62 </td>
   <td style="text-align:left;"> Median : 19.2 </td>
   <td style="text-align:left;"> Median :  8.8 </td>
   <td style="text-align:left;"> Median : 13.7 </td>
   <td style="text-align:left;"> Median : 21.4 </td>
   <td style="text-align:left;"> Median :1.47 </td>
   <td style="text-align:left;"> Median : 19 </td>
   <td style="text-align:left;"> Median : 15.6 </td>
   <td style="text-align:left;"> Median :1.1e+01 </td>
   <td style="text-align:left;"> Median :   10 </td>
   <td style="text-align:left;"> Median :   10 </td>
   <td style="text-align:left;"> Median :   8.28 </td>
   <td style="text-align:left;"> Median :  5.415 </td>
   <td style="text-align:left;"> Median :0.631 </td>
   <td style="text-align:left;"> Median :0.10491 </td>
   <td style="text-align:left;"> Median : 0.94 </td>
   <td style="text-align:left;"> Median :  7.6 </td>
   <td style="text-align:left;"> Median :0.09301 </td>
   <td style="text-align:left;"> Median : 1.0 </td>
   <td style="text-align:left;"> Median :3.9e+04 </td>
   <td style="text-align:left;"> Median :3.8e+08 </td>
   <td style="text-align:left;"> Median :1.0e+12 </td>
   <td style="text-align:left;"> Median : 3.1e+09 </td>
   <td style="text-align:left;"> Median : 0.33 </td>
   <td style="text-align:left;"> Median : 0.027 </td>
   <td style="text-align:left;"> Median :0.00061 </td>
   <td style="text-align:left;"> Median :3.0e-03 </td>
   <td style="text-align:left;"> Median :2.1e+06 </td>
   <td style="text-align:left;"> Median :6.8e+11 </td>
   <td style="text-align:left;"> Median :1.8e+12 </td>
   <td style="text-align:left;"> Median : 1.7e+09 </td>
   <td style="text-align:left;"> Median :3.4e+08 </td>
   <td style="text-align:left;"> Median :3.7e+10 </td>
   <td style="text-align:left;"> Median :1.7e-01 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :1.0e-03 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :  1035 </td>
   <td style="text-align:left;"> Median : 5.0e-05 </td>
   <td style="text-align:left;"> Median :2.0e-01 </td>
   <td style="text-align:left;"> Median :4.2e-01 </td>
   <td style="text-align:left;"> Median :0.62 </td>
   <td style="text-align:left;"> Median :  4.62 </td>
   <td style="text-align:left;"> Median : 2.89 </td>
   <td style="text-align:left;"> Median : 2.19 </td>
   <td style="text-align:left;"> Median :2.4e+04 </td>
   <td style="text-align:left;"> Median :1.4e+08 </td>
   <td style="text-align:left;"> Median :2.3e+11 </td>
   <td style="text-align:left;"> Median : 1.0e+09 </td>
   <td style="text-align:left;"> Median :  0.502 </td>
   <td style="text-align:left;"> Median :6.4e-02 </td>
   <td style="text-align:left;"> Median :2.2e-03 </td>
   <td style="text-align:left;"> Median :1.0e-02 </td>
   <td style="text-align:left;"> Median :1.3e+06 </td>
   <td style="text-align:left;"> Median :2.4e+11 </td>
   <td style="text-align:left;"> Median :6.0e+11 </td>
   <td style="text-align:left;"> Median : 5.8e+08 </td>
   <td style="text-align:left;"> Median :1.6e+08 </td>
   <td style="text-align:left;"> Median :1.4e+10 </td>
   <td style="text-align:left;"> Median :3.8e-01 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :0.0e+00 </td>
   <td style="text-align:left;"> Median :  676.9 </td>
   <td style="text-align:left;"> Median : 1.0e-05 </td>
   <td style="text-align:left;"> Median :2.1e-01 </td>
   <td style="text-align:left;"> Median :4.2e-01 </td>
   <td style="text-align:left;"> Median :0.6239 </td>
   <td style="text-align:left;"> Median :  4.33 </td>
   <td style="text-align:left;"> Median : 2.72 </td>
   <td style="text-align:left;"> Median : 2.08 </td>
   <td style="text-align:left;"> Median : 18.2 </td>
   <td style="text-align:left;"> Median : 15.7 </td>
   <td style="text-align:left;"> Median :  8.40 </td>
   <td style="text-align:left;"> Median : 12.0 </td>
   <td style="text-align:left;"> Median :  7.46 </td>
   <td style="text-align:left;"> Median : 17 </td>
   <td style="text-align:left;"> Median : 14.03 </td>
   <td style="text-align:left;"> Median : 13.4 </td>
   <td style="text-align:left;"> Median :6.9e+00 </td>
   <td style="text-align:left;"> Median : 20.8 </td>
   <td style="text-align:left;"> Median : 21.4 </td>
   <td style="text-align:left;"> Median : 15.390 </td>
   <td style="text-align:left;"> Median : 18.7 </td>
   <td style="text-align:left;"> Median :  7.76 </td>
   <td style="text-align:left;"> Median : 13.7 </td>
   <td style="text-align:left;"> Median : 21.3 </td>
   <td style="text-align:left;"> Median :1.59 </td>
   <td style="text-align:left;"> Median : 17.66 </td>
   <td style="text-align:left;"> Median : 15.6 </td>
   <td style="text-align:left;"> Median :9.4e+00 </td>
   <td style="text-align:left;"> Median : 14.2 </td>
   <td style="text-align:left;"> Median : 12.71 </td>
   <td style="text-align:left;"> Median :  8.0 </td>
   <td style="text-align:left;"> Median : 10.2 </td>
   <td style="text-align:left;"> Median : 6.11 </td>
   <td style="text-align:left;"> Median : 13.4 </td>
   <td style="text-align:left;"> Median : 11.06 </td>
   <td style="text-align:left;"> Median : 10.2 </td>
   <td style="text-align:left;"> Median :5.1e+00 </td>
   <td style="text-align:left;"> Median : 17.45 </td>
   <td style="text-align:left;"> Median : 16.6 </td>
   <td style="text-align:left;"> Median : 13.383 </td>
   <td style="text-align:left;"> Median : 14.85 </td>
   <td style="text-align:left;"> Median :  7.34 </td>
   <td style="text-align:left;"> Median : 11.1 </td>
   <td style="text-align:left;"> Median : 17.77 </td>
   <td style="text-align:left;"> Median :1.59 </td>
   <td style="text-align:left;"> Median : 15.4 </td>
   <td style="text-align:left;"> Median : 12.5 </td>
   <td style="text-align:left;"> Median :8.8e+00 </td>
   <td style="text-align:left;"> Median :2.05 </td>
   <td style="text-align:left;"> Median : 9.2e-13 </td>
   <td style="text-align:left;"> Median :0.124 </td>
   <td style="text-align:left;"> Median :0.0154 </td>
   <td style="text-align:left;"> Median :-0.675 </td>
   <td style="text-align:left;"> Median : 1.903 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Mean   : 13 </td>
   <td style="text-align:left;"> Mean   :  97 </td>
   <td style="text-align:left;"> Mean   : 13 </td>
   <td style="text-align:left;"> Mean   : 100 </td>
   <td style="text-align:left;"> Mean   :  897 </td>
   <td style="text-align:left;"> Mean   : 18.78 </td>
   <td style="text-align:left;"> Mean   :0.02477 </td>
   <td style="text-align:left;"> Mean   :0.1282 </td>
   <td style="text-align:left;"> Mean   :0 </td>
   <td style="text-align:left;"> Mean   : 1.42 </td>
   <td style="text-align:left;"> Mean   :0.2326 </td>
   <td style="text-align:left;"> Mean   :   355 </td>
   <td style="text-align:left;"> Mean   :   355 </td>
   <td style="text-align:left;"> Mean   :  34.68 </td>
   <td style="text-align:left;"> Mean   : 18.57 </td>
   <td style="text-align:left;"> Mean   :0.621 </td>
   <td style="text-align:left;"> Mean   :0.223 </td>
   <td style="text-align:left;"> Mean   : 1.42 </td>
   <td style="text-align:left;"> Mean   : 10.2 </td>
   <td style="text-align:left;"> Mean   :0.2298 </td>
   <td style="text-align:left;"> Mean   : 1.1 </td>
   <td style="text-align:left;"> Mean   :1.8e+06 </td>
   <td style="text-align:left;"> Mean   :1.2e+13 </td>
   <td style="text-align:left;"> Mean   :2.3e+20 </td>
   <td style="text-align:left;"> Mean   :4.8e+16 </td>
   <td style="text-align:left;"> Mean   : 0.49 </td>
   <td style="text-align:left;"> Mean   :0.062 </td>
   <td style="text-align:left;"> Mean   :0.00203 </td>
   <td style="text-align:left;"> Mean   :5.7e-02 </td>
   <td style="text-align:left;"> Mean   :3.6e+09 </td>
   <td style="text-align:left;"> Mean   :2.2e+20 </td>
   <td style="text-align:left;"> Mean   :1.6e+23 </td>
   <td style="text-align:left;"> Mean   :3.1e+16 </td>
   <td style="text-align:left;"> Mean   :1.9e+16 </td>
   <td style="text-align:left;"> Mean   :2.5e+18 </td>
   <td style="text-align:left;"> Mean   :5.5e-01 </td>
   <td style="text-align:left;"> Mean   :9.1e-02 </td>
   <td style="text-align:left;"> Mean   :     22 </td>
   <td style="text-align:left;"> Mean   :3.9e-02 </td>
   <td style="text-align:left;"> Mean   :2.6e-02 </td>
   <td style="text-align:left;"> Mean   :7.4e-01 </td>
   <td style="text-align:left;"> Mean   :  4335 </td>
   <td style="text-align:left;"> Mean   : 4.7e-02 </td>
   <td style="text-align:left;"> Mean   :0.24863 </td>
   <td style="text-align:left;"> Mean   :0.43057 </td>
   <td style="text-align:left;"> Mean   :0.557 </td>
   <td style="text-align:left;"> Mean   :  8.2 </td>
   <td style="text-align:left;"> Mean   : 4.6 </td>
   <td style="text-align:left;"> Mean   : 3.03 </td>
   <td style="text-align:left;"> Mean   :8.5e+05 </td>
   <td style="text-align:left;"> Mean   :1.5e+12 </td>
   <td style="text-align:left;"> Mean   :1.7e+18 </td>
   <td style="text-align:left;"> Mean   :2.7e+15 </td>
   <td style="text-align:left;"> Mean   :  0.730 </td>
   <td style="text-align:left;"> Mean   :0.14206 </td>
   <td style="text-align:left;"> Mean   :7.4e-03 </td>
   <td style="text-align:left;"> Mean   :2.8e-01 </td>
   <td style="text-align:left;"> Mean   :1.1e+09 </td>
   <td style="text-align:left;"> Mean   :8.4e+18 </td>
   <td style="text-align:left;"> Mean   :1.1e+21 </td>
   <td style="text-align:left;"> Mean   :1.7e+15 </td>
   <td style="text-align:left;"> Mean   :1.0e+15 </td>
   <td style="text-align:left;"> Mean   :3.2e+16 </td>
   <td style="text-align:left;"> Mean   :1.8e+00 </td>
   <td style="text-align:left;"> Mean   :6.9e-01 </td>
   <td style="text-align:left;"> Mean   :2.1e+04 </td>
   <td style="text-align:left;"> Mean   :2.1e-01 </td>
   <td style="text-align:left;"> Mean   :1.7e-01 </td>
   <td style="text-align:left;"> Mean   :7.3e+01 </td>
   <td style="text-align:left;"> Mean   : 2321 </td>
   <td style="text-align:left;"> Mean   : 4.8e-02 </td>
   <td style="text-align:left;"> Mean   :2.5e-01 </td>
   <td style="text-align:left;"> Mean   :0.43089 </td>
   <td style="text-align:left;"> Mean   :0.559 </td>
   <td style="text-align:left;"> Mean   :  7.9 </td>
   <td style="text-align:left;"> Mean   : 4.3 </td>
   <td style="text-align:left;"> Mean   : 2.86 </td>
   <td style="text-align:left;"> Mean   : 42.4 </td>
   <td style="text-align:left;"> Mean   : 27.1 </td>
   <td style="text-align:left;"> Mean   : 17.89 </td>
   <td style="text-align:left;"> Mean   : 29.1 </td>
   <td style="text-align:left;"> Mean   : 15.56 </td>
   <td style="text-align:left;"> Mean   : 36.2 </td>
   <td style="text-align:left;"> Mean   : 32.7 </td>
   <td style="text-align:left;"> Mean   : 25.2 </td>
   <td style="text-align:left;"> Mean   : 15.317 </td>
   <td style="text-align:left;"> Mean   : 39.6 </td>
   <td style="text-align:left;"> Mean   : 48.1 </td>
   <td style="text-align:left;"> Mean   : 29.1 </td>
   <td style="text-align:left;"> Mean   : 43 </td>
   <td style="text-align:left;"> Mean   : 18.88 </td>
   <td style="text-align:left;"> Mean   : 29.9 </td>
   <td style="text-align:left;"> Mean   : 45.2 </td>
   <td style="text-align:left;"> Mean   :1.42 </td>
   <td style="text-align:left;"> Mean   : 38.9 </td>
   <td style="text-align:left;"> Mean   : 37.7 </td>
   <td style="text-align:left;"> Mean   : 21.21 </td>
   <td style="text-align:left;"> Mean   : 31.5 </td>
   <td style="text-align:left;"> Mean   : 19.9 </td>
   <td style="text-align:left;"> Mean   : 15.38 </td>
   <td style="text-align:left;"> Mean   : 23.2 </td>
   <td style="text-align:left;"> Mean   :11.73 </td>
   <td style="text-align:left;"> Mean   : 26.6 </td>
   <td style="text-align:left;"> Mean   : 25.8 </td>
   <td style="text-align:left;"> Mean   : 18.0 </td>
   <td style="text-align:left;"> Mean   :1.3e+01 </td>
   <td style="text-align:left;"> Mean   : 29.2 </td>
   <td style="text-align:left;"> Mean   : 35.6 </td>
   <td style="text-align:left;"> Mean   : 22.6 </td>
   <td style="text-align:left;"> Mean   : 32.8 </td>
   <td style="text-align:left;"> Mean   : 15.47 </td>
   <td style="text-align:left;"> Mean   : 22.8 </td>
   <td style="text-align:left;"> Mean   : 34 </td>
   <td style="text-align:left;"> Mean   :1.41 </td>
   <td style="text-align:left;"> Mean   : 29.74 </td>
   <td style="text-align:left;"> Mean   : 28.7 </td>
   <td style="text-align:left;"> Mean   :1.8e+01 </td>
   <td style="text-align:left;"> Mean   :  298 </td>
   <td style="text-align:left;"> Mean   :  298 </td>
   <td style="text-align:left;"> Mean   :  26.83 </td>
   <td style="text-align:left;"> Mean   : 16.035 </td>
   <td style="text-align:left;"> Mean   :0.677 </td>
   <td style="text-align:left;"> Mean   :0.2128 </td>
   <td style="text-align:left;"> Mean   : 1.42 </td>
   <td style="text-align:left;"> Mean   : 10.2 </td>
   <td style="text-align:left;"> Mean   :0.2166 </td>
   <td style="text-align:left;"> Mean   : 1.3 </td>
   <td style="text-align:left;"> Mean   :1.4e+06 </td>
   <td style="text-align:left;"> Mean   :6.2e+12 </td>
   <td style="text-align:left;"> Mean   :6.4e+19 </td>
   <td style="text-align:left;"> Mean   :2.3e+16 </td>
   <td style="text-align:left;"> Mean   : 0.53 </td>
   <td style="text-align:left;"> Mean   : 0.073 </td>
   <td style="text-align:left;"> Mean   :0.00262 </td>
   <td style="text-align:left;"> Mean   :1.3e-01 </td>
   <td style="text-align:left;"> Mean   :2.7e+09 </td>
   <td style="text-align:left;"> Mean   :9.3e+19 </td>
   <td style="text-align:left;"> Mean   :1.0e+23 </td>
   <td style="text-align:left;"> Mean   :1.5e+16 </td>
   <td style="text-align:left;"> Mean   :9.0e+15 </td>
   <td style="text-align:left;"> Mean   :1.6e+18 </td>
   <td style="text-align:left;"> Mean   :7.4e-01 </td>
   <td style="text-align:left;"> Mean   :2.1e-01 </td>
   <td style="text-align:left;"> Mean   :     46 </td>
   <td style="text-align:left;"> Mean   :1.0e-01 </td>
   <td style="text-align:left;"> Mean   :  0.084 </td>
   <td style="text-align:left;"> Mean   :1.3e+00 </td>
   <td style="text-align:left;"> Mean   :  3354 </td>
   <td style="text-align:left;"> Mean   : 4.1e-02 </td>
   <td style="text-align:left;"> Mean   :2.6e-01 </td>
   <td style="text-align:left;"> Mean   :0.43148 </td>
   <td style="text-align:left;"> Mean   :0.566 </td>
   <td style="text-align:left;"> Mean   :  7.7 </td>
   <td style="text-align:left;"> Mean   : 4.16 </td>
   <td style="text-align:left;"> Mean   : 2.76 </td>
   <td style="text-align:left;"> Mean   :7.1e+05 </td>
   <td style="text-align:left;"> Mean   :1.1e+12 </td>
   <td style="text-align:left;"> Mean   :9.6e+17 </td>
   <td style="text-align:left;"> Mean   :1.9e+15 </td>
   <td style="text-align:left;"> Mean   :  0.753 </td>
   <td style="text-align:left;"> Mean   :1.5e-01 </td>
   <td style="text-align:left;"> Mean   :8.1e-03 </td>
   <td style="text-align:left;"> Mean   :6.1e-01 </td>
   <td style="text-align:left;"> Mean   :9.4e+08 </td>
   <td style="text-align:left;"> Mean   :5.6e+18 </td>
   <td style="text-align:left;"> Mean   :8.3e+20 </td>
   <td style="text-align:left;"> Mean   :1.2e+15 </td>
   <td style="text-align:left;"> Mean   :7.7e+14 </td>
   <td style="text-align:left;"> Mean   :2.3e+16 </td>
   <td style="text-align:left;"> Mean   :2.3e+00 </td>
   <td style="text-align:left;"> Mean   :1.6e+00 </td>
   <td style="text-align:left;"> Mean   :3.8e+04 </td>
   <td style="text-align:left;"> Mean   :5.4e-01 </td>
   <td style="text-align:left;"> Mean   :4.8e-01 </td>
   <td style="text-align:left;"> Mean   :1.2e+02 </td>
   <td style="text-align:left;"> Mean   : 2004 </td>
   <td style="text-align:left;"> Mean   : 4.2e-02 </td>
   <td style="text-align:left;"> Mean   :2.6e-01 </td>
   <td style="text-align:left;"> Mean   :0.43219 </td>
   <td style="text-align:left;"> Mean   :0.568 </td>
   <td style="text-align:left;"> Mean   :  7.4 </td>
   <td style="text-align:left;"> Mean   : 3.97 </td>
   <td style="text-align:left;"> Mean   : 2.62 </td>
   <td style="text-align:left;"> Mean   : 37.0 </td>
   <td style="text-align:left;"> Mean   : 23 </td>
   <td style="text-align:left;"> Mean   : 16.33 </td>
   <td style="text-align:left;"> Mean   : 25.7 </td>
   <td style="text-align:left;"> Mean   : 13.83 </td>
   <td style="text-align:left;"> Mean   : 31.3 </td>
   <td style="text-align:left;"> Mean   : 28.4 </td>
   <td style="text-align:left;"> Mean   : 22.1 </td>
   <td style="text-align:left;"> Mean   : 13.608 </td>
   <td style="text-align:left;"> Mean   : 33.9 </td>
   <td style="text-align:left;"> Mean   : 41.6 </td>
   <td style="text-align:left;"> Mean   : 24.79 </td>
   <td style="text-align:left;"> Mean   : 37.4 </td>
   <td style="text-align:left;"> Mean   : 16.85 </td>
   <td style="text-align:left;"> Mean   : 25.8 </td>
   <td style="text-align:left;"> Mean   : 38.8 </td>
   <td style="text-align:left;"> Mean   :1.5 </td>
   <td style="text-align:left;"> Mean   : 33 </td>
   <td style="text-align:left;"> Mean   : 33 </td>
   <td style="text-align:left;"> Mean   : 18.103 </td>
   <td style="text-align:left;"> Mean   : 28.9 </td>
   <td style="text-align:left;"> Mean   : 18.1 </td>
   <td style="text-align:left;"> Mean   : 14.6 </td>
   <td style="text-align:left;"> Mean   : 21.4 </td>
   <td style="text-align:left;"> Mean   :10.99 </td>
   <td style="text-align:left;"> Mean   : 24.3 </td>
   <td style="text-align:left;"> Mean   : 23.05 </td>
   <td style="text-align:left;"> Mean   : 16.6 </td>
   <td style="text-align:left;"> Mean   : 11.866 </td>
   <td style="text-align:left;"> Mean   : 26.4 </td>
   <td style="text-align:left;"> Mean   : 32.2 </td>
   <td style="text-align:left;"> Mean   : 20.21 </td>
   <td style="text-align:left;"> Mean   : 29.4 </td>
   <td style="text-align:left;"> Mean   : 14.4 </td>
   <td style="text-align:left;"> Mean   : 20.7 </td>
   <td style="text-align:left;"> Mean   : 30.1 </td>
   <td style="text-align:left;"> Mean   :1.50 </td>
   <td style="text-align:left;"> Mean   : 26 </td>
   <td style="text-align:left;"> Mean   : 26.2 </td>
   <td style="text-align:left;"> Mean   :1.6e+01 </td>
   <td style="text-align:left;"> Mean   :  249 </td>
   <td style="text-align:left;"> Mean   :  249 </td>
   <td style="text-align:left;"> Mean   :  20.75 </td>
   <td style="text-align:left;"> Mean   : 13.723 </td>
   <td style="text-align:left;"> Mean   :0.731 </td>
   <td style="text-align:left;"> Mean   :0.20179 </td>
   <td style="text-align:left;"> Mean   : 1.42 </td>
   <td style="text-align:left;"> Mean   : 10.2 </td>
   <td style="text-align:left;"> Mean   :0.20337 </td>
   <td style="text-align:left;"> Mean   : 1.4 </td>
   <td style="text-align:left;"> Mean   :1.0e+06 </td>
   <td style="text-align:left;"> Mean   :3.3e+12 </td>
   <td style="text-align:left;"> Mean   :2.0e+19 </td>
   <td style="text-align:left;"> Mean   : 1.1e+16 </td>
   <td style="text-align:left;"> Mean   : 0.57 </td>
   <td style="text-align:left;"> Mean   : 0.087 </td>
   <td style="text-align:left;"> Mean   :0.00349 </td>
   <td style="text-align:left;"> Mean   :3.7e-01 </td>
   <td style="text-align:left;"> Mean   :2.0e+09 </td>
   <td style="text-align:left;"> Mean   :4.3e+19 </td>
   <td style="text-align:left;"> Mean   :6.4e+22 </td>
   <td style="text-align:left;"> Mean   : 7.2e+15 </td>
   <td style="text-align:left;"> Mean   :4.5e+15 </td>
   <td style="text-align:left;"> Mean   :9.7e+17 </td>
   <td style="text-align:left;"> Mean   :1.1e+00 </td>
   <td style="text-align:left;"> Mean   :9.7e-01 </td>
   <td style="text-align:left;"> Mean   :3.5e+02 </td>
   <td style="text-align:left;"> Mean   :3.4e-01 </td>
   <td style="text-align:left;"> Mean   :3.1e-01 </td>
   <td style="text-align:left;"> Mean   :4.6e+00 </td>
   <td style="text-align:left;"> Mean   :  2594 </td>
   <td style="text-align:left;"> Mean   : 3.6e-02 </td>
   <td style="text-align:left;"> Mean   :2.7e-01 </td>
   <td style="text-align:left;"> Mean   :4.4e-01 </td>
   <td style="text-align:left;"> Mean   :0.58 </td>
   <td style="text-align:left;"> Mean   :  7.03 </td>
   <td style="text-align:left;"> Mean   : 3.76 </td>
   <td style="text-align:left;"> Mean   : 2.49 </td>
   <td style="text-align:left;"> Mean   :5.9e+05 </td>
   <td style="text-align:left;"> Mean   :7.6e+11 </td>
   <td style="text-align:left;"> Mean   :5.7e+17 </td>
   <td style="text-align:left;"> Mean   : 1.4e+15 </td>
   <td style="text-align:left;"> Mean   :  0.771 </td>
   <td style="text-align:left;"> Mean   :1.6e-01 </td>
   <td style="text-align:left;"> Mean   :9.5e-03 </td>
   <td style="text-align:left;"> Mean   :2.8e+00 </td>
   <td style="text-align:left;"> Mean   :7.7e+08 </td>
   <td style="text-align:left;"> Mean   :3.8e+18 </td>
   <td style="text-align:left;"> Mean   :5.8e+20 </td>
   <td style="text-align:left;"> Mean   : 8.9e+14 </td>
   <td style="text-align:left;"> Mean   :5.7e+14 </td>
   <td style="text-align:left;"> Mean   :1.7e+16 </td>
   <td style="text-align:left;"> Mean   :3.8e+00 </td>
   <td style="text-align:left;"> Mean   :1.6e+01 </td>
   <td style="text-align:left;"> Mean   :3.0e+05 </td>
   <td style="text-align:left;"> Mean   :3.0e+00 </td>
   <td style="text-align:left;"> Mean   :3.1e+00 </td>
   <td style="text-align:left;"> Mean   :5.0e+02 </td>
   <td style="text-align:left;"> Mean   : 1715.4 </td>
   <td style="text-align:left;"> Mean   : 3.4e-02 </td>
   <td style="text-align:left;"> Mean   :2.7e-01 </td>
   <td style="text-align:left;"> Mean   :4.4e-01 </td>
   <td style="text-align:left;"> Mean   :0.5825 </td>
   <td style="text-align:left;"> Mean   :  6.80 </td>
   <td style="text-align:left;"> Mean   : 3.60 </td>
   <td style="text-align:left;"> Mean   : 2.38 </td>
   <td style="text-align:left;"> Mean   : 32.2 </td>
   <td style="text-align:left;"> Mean   : 19.9 </td>
   <td style="text-align:left;"> Mean   : 15.05 </td>
   <td style="text-align:left;"> Mean   : 22.8 </td>
   <td style="text-align:left;"> Mean   : 12.26 </td>
   <td style="text-align:left;"> Mean   : 27 </td>
   <td style="text-align:left;"> Mean   : 24.12 </td>
   <td style="text-align:left;"> Mean   : 19.1 </td>
   <td style="text-align:left;"> Mean   :1.2e+01 </td>
   <td style="text-align:left;"> Mean   : 28.7 </td>
   <td style="text-align:left;"> Mean   : 35.5 </td>
   <td style="text-align:left;"> Mean   : 20.821 </td>
   <td style="text-align:left;"> Mean   : 31.7 </td>
   <td style="text-align:left;"> Mean   : 15.09 </td>
   <td style="text-align:left;"> Mean   : 22.3 </td>
   <td style="text-align:left;"> Mean   : 32.7 </td>
   <td style="text-align:left;"> Mean   :1.66 </td>
   <td style="text-align:left;"> Mean   : 27.72 </td>
   <td style="text-align:left;"> Mean   : 28.5 </td>
   <td style="text-align:left;"> Mean   :1.5e+01 </td>
   <td style="text-align:left;"> Mean   : 26.6 </td>
   <td style="text-align:left;"> Mean   : 16.21 </td>
   <td style="text-align:left;"> Mean   : 14.0 </td>
   <td style="text-align:left;"> Mean   : 19.9 </td>
   <td style="text-align:left;"> Mean   :10.30 </td>
   <td style="text-align:left;"> Mean   : 22.1 </td>
   <td style="text-align:left;"> Mean   : 20.31 </td>
   <td style="text-align:left;"> Mean   : 15.3 </td>
   <td style="text-align:left;"> Mean   :1.1e+01 </td>
   <td style="text-align:left;"> Mean   : 23.50 </td>
   <td style="text-align:left;"> Mean   : 28.6 </td>
   <td style="text-align:left;"> Mean   : 17.850 </td>
   <td style="text-align:left;"> Mean   : 25.98 </td>
   <td style="text-align:left;"> Mean   : 13.49 </td>
   <td style="text-align:left;"> Mean   : 18.8 </td>
   <td style="text-align:left;"> Mean   : 26.57 </td>
   <td style="text-align:left;"> Mean   :1.65 </td>
   <td style="text-align:left;"> Mean   : 23.1 </td>
   <td style="text-align:left;"> Mean   : 24.0 </td>
   <td style="text-align:left;"> Mean   :1.4e+01 </td>
   <td style="text-align:left;"> Mean   :2.13 </td>
   <td style="text-align:left;"> Mean   : 3.3e-11 </td>
   <td style="text-align:left;"> Mean   :0.130 </td>
   <td style="text-align:left;"> Mean   :0.0198 </td>
   <td style="text-align:left;"> Mean   :-0.708 </td>
   <td style="text-align:left;"> Mean   : 2.678 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> 3rd Qu.: 17 </td>
   <td style="text-align:left;"> 3rd Qu.: 121 </td>
   <td style="text-align:left;"> 3rd Qu.: 18 </td>
   <td style="text-align:left;"> 3rd Qu.: 128 </td>
   <td style="text-align:left;"> 3rd Qu.:  835 </td>
   <td style="text-align:left;"> 3rd Qu.: 20.89 </td>
   <td style="text-align:left;"> 3rd Qu.:0.02965 </td>
   <td style="text-align:left;"> 3rd Qu.:0.1480 </td>
   <td style="text-align:left;"> 3rd Qu.:0 </td>
   <td style="text-align:left;"> 3rd Qu.: 1.58 </td>
   <td style="text-align:left;"> 3rd Qu.:0.2650 </td>
   <td style="text-align:left;"> 3rd Qu.:   155 </td>
   <td style="text-align:left;"> 3rd Qu.:   155 </td>
   <td style="text-align:left;"> 3rd Qu.:  36.66 </td>
   <td style="text-align:left;"> 3rd Qu.: 20.23 </td>
   <td style="text-align:left;"> 3rd Qu.:0.737 </td>
   <td style="text-align:left;"> 3rd Qu.:0.247 </td>
   <td style="text-align:left;"> 3rd Qu.: 1.58 </td>
   <td style="text-align:left;"> 3rd Qu.: 11.8 </td>
   <td style="text-align:left;"> 3rd Qu.:0.2498 </td>
   <td style="text-align:left;"> 3rd Qu.: 1.0 </td>
   <td style="text-align:left;"> 3rd Qu.:6.7e+05 </td>
   <td style="text-align:left;"> 3rd Qu.:9.3e+10 </td>
   <td style="text-align:left;"> 3rd Qu.:3.0e+15 </td>
   <td style="text-align:left;"> 3rd Qu.:7.1e+12 </td>
   <td style="text-align:left;"> 3rd Qu.: 0.61 </td>
   <td style="text-align:left;"> 3rd Qu.:0.070 </td>
   <td style="text-align:left;"> 3rd Qu.:0.00182 </td>
   <td style="text-align:left;"> 3rd Qu.:3.2e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:1.6e+08 </td>
   <td style="text-align:left;"> 3rd Qu.:2.8e+15 </td>
   <td style="text-align:left;"> 3rd Qu.:1.2e+16 </td>
   <td style="text-align:left;"> 3rd Qu.:4.3e+12 </td>
   <td style="text-align:left;"> 3rd Qu.:1.8e+12 </td>
   <td style="text-align:left;"> 3rd Qu.:5.4e+13 </td>
   <td style="text-align:left;"> 3rd Qu.:5.9e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:3.6e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:      0 </td>
   <td style="text-align:left;"> 3rd Qu.:2.0e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:9.6e-03 </td>
   <td style="text-align:left;"> 3rd Qu.:2.0e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:  4583 </td>
   <td style="text-align:left;"> 3rd Qu.: 8.3e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:0.36817 </td>
   <td style="text-align:left;"> 3rd Qu.:0.62528 </td>
   <td style="text-align:left;"> 3rd Qu.:0.748 </td>
   <td style="text-align:left;"> 3rd Qu.: 10.3 </td>
   <td style="text-align:left;"> 3rd Qu.: 5.4 </td>
   <td style="text-align:left;"> 3rd Qu.: 3.54 </td>
   <td style="text-align:left;"> 3rd Qu.:3.0e+05 </td>
   <td style="text-align:left;"> 3rd Qu.:1.8e+10 </td>
   <td style="text-align:left;"> 3rd Qu.:2.5e+14 </td>
   <td style="text-align:left;"> 3rd Qu.:1.5e+12 </td>
   <td style="text-align:left;"> 3rd Qu.:  0.940 </td>
   <td style="text-align:left;"> 3rd Qu.:0.17312 </td>
   <td style="text-align:left;"> 3rd Qu.:7.3e-03 </td>
   <td style="text-align:left;"> 3rd Qu.:1.2e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:6.8e+07 </td>
   <td style="text-align:left;"> 3rd Qu.:5.1e+14 </td>
   <td style="text-align:left;"> 3rd Qu.:2.4e+15 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e+12 </td>
   <td style="text-align:left;"> 3rd Qu.:5.3e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:1.1e+13 </td>
   <td style="text-align:left;"> 3rd Qu.:1.5e+00 </td>
   <td style="text-align:left;"> 3rd Qu.:2.2e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e+00 </td>
   <td style="text-align:left;"> 3rd Qu.:8.3e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:4.7e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e+00 </td>
   <td style="text-align:left;"> 3rd Qu.: 2529 </td>
   <td style="text-align:left;"> 3rd Qu.: 9.0e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:3.8e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:0.63050 </td>
   <td style="text-align:left;"> 3rd Qu.:0.754 </td>
   <td style="text-align:left;"> 3rd Qu.: 10.0 </td>
   <td style="text-align:left;"> 3rd Qu.: 5.1 </td>
   <td style="text-align:left;"> 3rd Qu.: 3.30 </td>
   <td style="text-align:left;"> 3rd Qu.: 55.2 </td>
   <td style="text-align:left;"> 3rd Qu.: 33.1 </td>
   <td style="text-align:left;"> 3rd Qu.: 22.93 </td>
   <td style="text-align:left;"> 3rd Qu.: 38.2 </td>
   <td style="text-align:left;"> 3rd Qu.: 20.23 </td>
   <td style="text-align:left;"> 3rd Qu.: 46.6 </td>
   <td style="text-align:left;"> 3rd Qu.: 43.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 31.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 19.967 </td>
   <td style="text-align:left;"> 3rd Qu.: 49.5 </td>
   <td style="text-align:left;"> 3rd Qu.: 63.2 </td>
   <td style="text-align:left;"> 3rd Qu.: 36.5 </td>
   <td style="text-align:left;"> 3rd Qu.: 58 </td>
   <td style="text-align:left;"> 3rd Qu.: 25.26 </td>
   <td style="text-align:left;"> 3rd Qu.: 38.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 58.7 </td>
   <td style="text-align:left;"> 3rd Qu.:1.56 </td>
   <td style="text-align:left;"> 3rd Qu.: 51.1 </td>
   <td style="text-align:left;"> 3rd Qu.: 49.4 </td>
   <td style="text-align:left;"> 3rd Qu.: 28.20 </td>
   <td style="text-align:left;"> 3rd Qu.: 40.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 24.6 </td>
   <td style="text-align:left;"> 3rd Qu.: 19.73 </td>
   <td style="text-align:left;"> 3rd Qu.: 30.3 </td>
   <td style="text-align:left;"> 3rd Qu.:15.02 </td>
   <td style="text-align:left;"> 3rd Qu.: 34.4 </td>
   <td style="text-align:left;"> 3rd Qu.: 34.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 22.7 </td>
   <td style="text-align:left;"> 3rd Qu.:1.8e+01 </td>
   <td style="text-align:left;"> 3rd Qu.: 36.0 </td>
   <td style="text-align:left;"> 3rd Qu.: 46.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 28.5 </td>
   <td style="text-align:left;"> 3rd Qu.: 43.6 </td>
   <td style="text-align:left;"> 3rd Qu.: 20.29 </td>
   <td style="text-align:left;"> 3rd Qu.: 29.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 43 </td>
   <td style="text-align:left;"> 3rd Qu.:1.55 </td>
   <td style="text-align:left;"> 3rd Qu.: 38.40 </td>
   <td style="text-align:left;"> 3rd Qu.: 37.5 </td>
   <td style="text-align:left;"> 3rd Qu.:2.4e+01 </td>
   <td style="text-align:left;"> 3rd Qu.:  106 </td>
   <td style="text-align:left;"> 3rd Qu.:  106 </td>
   <td style="text-align:left;"> 3rd Qu.:  27.44 </td>
   <td style="text-align:left;"> 3rd Qu.: 17.681 </td>
   <td style="text-align:left;"> 3rd Qu.:0.804 </td>
   <td style="text-align:left;"> 3rd Qu.:0.2358 </td>
   <td style="text-align:left;"> 3rd Qu.: 1.58 </td>
   <td style="text-align:left;"> 3rd Qu.: 11.8 </td>
   <td style="text-align:left;"> 3rd Qu.:0.2347 </td>
   <td style="text-align:left;"> 3rd Qu.: 1.0 </td>
   <td style="text-align:left;"> 3rd Qu.:4.3e+05 </td>
   <td style="text-align:left;"> 3rd Qu.:3.7e+10 </td>
   <td style="text-align:left;"> 3rd Qu.:7.5e+14 </td>
   <td style="text-align:left;"> 3rd Qu.:2.7e+12 </td>
   <td style="text-align:left;"> 3rd Qu.: 0.67 </td>
   <td style="text-align:left;"> 3rd Qu.: 0.080 </td>
   <td style="text-align:left;"> 3rd Qu.:0.00208 </td>
   <td style="text-align:left;"> 3rd Qu.:4.7e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:8.9e+07 </td>
   <td style="text-align:left;"> 3rd Qu.:8.6e+14 </td>
   <td style="text-align:left;"> 3rd Qu.:4.1e+15 </td>
   <td style="text-align:left;"> 3rd Qu.:1.8e+12 </td>
   <td style="text-align:left;"> 3rd Qu.:8.3e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:2.0e+13 </td>
   <td style="text-align:left;"> 3rd Qu.:7.5e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:5.0e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:      0 </td>
   <td style="text-align:left;"> 3rd Qu.:3.1e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:  0.016 </td>
   <td style="text-align:left;"> 3rd Qu.:2.8e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:  3430 </td>
   <td style="text-align:left;"> 3rd Qu.: 6.7e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:4.0e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:0.63990 </td>
   <td style="text-align:left;"> 3rd Qu.:0.762 </td>
   <td style="text-align:left;"> 3rd Qu.:  9.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 4.94 </td>
   <td style="text-align:left;"> 3rd Qu.: 3.28 </td>
   <td style="text-align:left;"> 3rd Qu.:2.1e+05 </td>
   <td style="text-align:left;"> 3rd Qu.:9.2e+09 </td>
   <td style="text-align:left;"> 3rd Qu.:9.3e+13 </td>
   <td style="text-align:left;"> 3rd Qu.:7.5e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:  0.962 </td>
   <td style="text-align:left;"> 3rd Qu.:1.7e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:7.0e-03 </td>
   <td style="text-align:left;"> 3rd Qu.:1.4e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:4.4e+07 </td>
   <td style="text-align:left;"> 3rd Qu.:2.0e+14 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e+15 </td>
   <td style="text-align:left;"> 3rd Qu.:5.2e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:2.9e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:5.0e+12 </td>
   <td style="text-align:left;"> 3rd Qu.:1.6e+00 </td>
   <td style="text-align:left;"> 3rd Qu.:2.3e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e+00 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:6.1e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e+00 </td>
   <td style="text-align:left;"> 3rd Qu.: 2210 </td>
   <td style="text-align:left;"> 3rd Qu.: 7.0e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:4.1e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:0.64446 </td>
   <td style="text-align:left;"> 3rd Qu.:0.767 </td>
   <td style="text-align:left;"> 3rd Qu.:  9.4 </td>
   <td style="text-align:left;"> 3rd Qu.: 4.68 </td>
   <td style="text-align:left;"> 3rd Qu.: 3.07 </td>
   <td style="text-align:left;"> 3rd Qu.: 48.1 </td>
   <td style="text-align:left;"> 3rd Qu.: 29 </td>
   <td style="text-align:left;"> 3rd Qu.: 20.38 </td>
   <td style="text-align:left;"> 3rd Qu.: 33.5 </td>
   <td style="text-align:left;"> 3rd Qu.: 18.06 </td>
   <td style="text-align:left;"> 3rd Qu.: 40.6 </td>
   <td style="text-align:left;"> 3rd Qu.: 37.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 27.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 17.616 </td>
   <td style="text-align:left;"> 3rd Qu.: 42.6 </td>
   <td style="text-align:left;"> 3rd Qu.: 54.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 31.34 </td>
   <td style="text-align:left;"> 3rd Qu.: 49.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 22.40 </td>
   <td style="text-align:left;"> 3rd Qu.: 33.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 50.5 </td>
   <td style="text-align:left;"> 3rd Qu.:1.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 44 </td>
   <td style="text-align:left;"> 3rd Qu.: 43 </td>
   <td style="text-align:left;"> 3rd Qu.: 24.236 </td>
   <td style="text-align:left;"> 3rd Qu.: 37.5 </td>
   <td style="text-align:left;"> 3rd Qu.: 23.0 </td>
   <td style="text-align:left;"> 3rd Qu.: 18.1 </td>
   <td style="text-align:left;"> 3rd Qu.: 27.7 </td>
   <td style="text-align:left;"> 3rd Qu.:14.11 </td>
   <td style="text-align:left;"> 3rd Qu.: 31.6 </td>
   <td style="text-align:left;"> 3rd Qu.: 31.76 </td>
   <td style="text-align:left;"> 3rd Qu.: 21.2 </td>
   <td style="text-align:left;"> 3rd Qu.: 16.326 </td>
   <td style="text-align:left;"> 3rd Qu.: 33.0 </td>
   <td style="text-align:left;"> 3rd Qu.: 43.0 </td>
   <td style="text-align:left;"> 3rd Qu.: 26.23 </td>
   <td style="text-align:left;"> 3rd Qu.: 39.6 </td>
   <td style="text-align:left;"> 3rd Qu.: 18.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 27.2 </td>
   <td style="text-align:left;"> 3rd Qu.: 39.3 </td>
   <td style="text-align:left;"> 3rd Qu.:1.69 </td>
   <td style="text-align:left;"> 3rd Qu.: 35 </td>
   <td style="text-align:left;"> 3rd Qu.: 34.1 </td>
   <td style="text-align:left;"> 3rd Qu.:2.2e+01 </td>
   <td style="text-align:left;"> 3rd Qu.:   71 </td>
   <td style="text-align:left;"> 3rd Qu.:   71 </td>
   <td style="text-align:left;"> 3rd Qu.:  20.67 </td>
   <td style="text-align:left;"> 3rd Qu.: 15.309 </td>
   <td style="text-align:left;"> 3rd Qu.:0.870 </td>
   <td style="text-align:left;"> 3rd Qu.:0.22478 </td>
   <td style="text-align:left;"> 3rd Qu.: 1.58 </td>
   <td style="text-align:left;"> 3rd Qu.: 11.8 </td>
   <td style="text-align:left;"> 3rd Qu.:0.21934 </td>
   <td style="text-align:left;"> 3rd Qu.: 1.0 </td>
   <td style="text-align:left;"> 3rd Qu.:2.6e+05 </td>
   <td style="text-align:left;"> 3rd Qu.:1.4e+10 </td>
   <td style="text-align:left;"> 3rd Qu.:1.8e+14 </td>
   <td style="text-align:left;"> 3rd Qu.: 9.7e+11 </td>
   <td style="text-align:left;"> 3rd Qu.: 0.71 </td>
   <td style="text-align:left;"> 3rd Qu.: 0.086 </td>
   <td style="text-align:left;"> 3rd Qu.:0.00222 </td>
   <td style="text-align:left;"> 3rd Qu.:5.7e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:4.8e+07 </td>
   <td style="text-align:left;"> 3rd Qu.:2.3e+14 </td>
   <td style="text-align:left;"> 3rd Qu.:1.2e+15 </td>
   <td style="text-align:left;"> 3rd Qu.: 6.4e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:3.1e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:6.6e+12 </td>
   <td style="text-align:left;"> 3rd Qu.:8.5e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:6.0e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:0.0e+00 </td>
   <td style="text-align:left;"> 3rd Qu.:3.9e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:2.1e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:3.0e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:  2584 </td>
   <td style="text-align:left;"> 3rd Qu.: 4.7e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:4.3e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:6.6e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:0.78 </td>
   <td style="text-align:left;"> 3rd Qu.:  9.05 </td>
   <td style="text-align:left;"> 3rd Qu.: 4.48 </td>
   <td style="text-align:left;"> 3rd Qu.: 3.04 </td>
   <td style="text-align:left;"> 3rd Qu.:1.5e+05 </td>
   <td style="text-align:left;"> 3rd Qu.:4.5e+09 </td>
   <td style="text-align:left;"> 3rd Qu.:3.3e+13 </td>
   <td style="text-align:left;"> 3rd Qu.: 3.4e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:  0.948 </td>
   <td style="text-align:left;"> 3rd Qu.:1.6e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:6.3e-03 </td>
   <td style="text-align:left;"> 3rd Qu.:1.4e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:2.7e+07 </td>
   <td style="text-align:left;"> 3rd Qu.:7.1e+13 </td>
   <td style="text-align:left;"> 3rd Qu.:3.9e+14 </td>
   <td style="text-align:left;"> 3rd Qu.: 2.4e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:1.4e+11 </td>
   <td style="text-align:left;"> 3rd Qu.:2.1e+12 </td>
   <td style="text-align:left;"> 3rd Qu.:1.6e+00 </td>
   <td style="text-align:left;"> 3rd Qu.:2.0e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e+00 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:6.0e-02 </td>
   <td style="text-align:left;"> 3rd Qu.:1.0e+00 </td>
   <td style="text-align:left;"> 3rd Qu.: 1913.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 4.9e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:4.4e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:6.6e-01 </td>
   <td style="text-align:left;"> 3rd Qu.:0.7801 </td>
   <td style="text-align:left;"> 3rd Qu.:  8.81 </td>
   <td style="text-align:left;"> 3rd Qu.: 4.25 </td>
   <td style="text-align:left;"> 3rd Qu.: 2.85 </td>
   <td style="text-align:left;"> 3rd Qu.: 41.3 </td>
   <td style="text-align:left;"> 3rd Qu.: 24.8 </td>
   <td style="text-align:left;"> 3rd Qu.: 17.66 </td>
   <td style="text-align:left;"> 3rd Qu.: 28.8 </td>
   <td style="text-align:left;"> 3rd Qu.: 15.90 </td>
   <td style="text-align:left;"> 3rd Qu.: 35 </td>
   <td style="text-align:left;"> 3rd Qu.: 31.94 </td>
   <td style="text-align:left;"> 3rd Qu.: 24.3 </td>
   <td style="text-align:left;"> 3rd Qu.:1.5e+01 </td>
   <td style="text-align:left;"> 3rd Qu.: 36.4 </td>
   <td style="text-align:left;"> 3rd Qu.: 46.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 26.736 </td>
   <td style="text-align:left;"> 3rd Qu.: 42.1 </td>
   <td style="text-align:left;"> 3rd Qu.: 19.47 </td>
   <td style="text-align:left;"> 3rd Qu.: 28.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 43.1 </td>
   <td style="text-align:left;"> 3rd Qu.:1.88 </td>
   <td style="text-align:left;"> 3rd Qu.: 36.91 </td>
   <td style="text-align:left;"> 3rd Qu.: 36.5 </td>
   <td style="text-align:left;"> 3rd Qu.:2.0e+01 </td>
   <td style="text-align:left;"> 3rd Qu.: 33.8 </td>
   <td style="text-align:left;"> 3rd Qu.: 21.37 </td>
   <td style="text-align:left;"> 3rd Qu.: 16.2 </td>
   <td style="text-align:left;"> 3rd Qu.: 24.8 </td>
   <td style="text-align:left;"> 3rd Qu.:13.03 </td>
   <td style="text-align:left;"> 3rd Qu.: 28.7 </td>
   <td style="text-align:left;"> 3rd Qu.: 28.47 </td>
   <td style="text-align:left;"> 3rd Qu.: 19.5 </td>
   <td style="text-align:left;"> 3rd Qu.:1.4e+01 </td>
   <td style="text-align:left;"> 3rd Qu.: 30.10 </td>
   <td style="text-align:left;"> 3rd Qu.: 38.9 </td>
   <td style="text-align:left;"> 3rd Qu.: 24.036 </td>
   <td style="text-align:left;"> 3rd Qu.: 35.62 </td>
   <td style="text-align:left;"> 3rd Qu.: 16.97 </td>
   <td style="text-align:left;"> 3rd Qu.: 24.5 </td>
   <td style="text-align:left;"> 3rd Qu.: 35.61 </td>
   <td style="text-align:left;"> 3rd Qu.:1.88 </td>
   <td style="text-align:left;"> 3rd Qu.: 31.3 </td>
   <td style="text-align:left;"> 3rd Qu.: 30.5 </td>
   <td style="text-align:left;"> 3rd Qu.:1.9e+01 </td>
   <td style="text-align:left;"> 3rd Qu.:2.40 </td>
   <td style="text-align:left;"> 3rd Qu.: 5.1e-11 </td>
   <td style="text-align:left;"> 3rd Qu.:0.161 </td>
   <td style="text-align:left;"> 3rd Qu.:0.0258 </td>
   <td style="text-align:left;"> 3rd Qu.:-0.508 </td>
   <td style="text-align:left;"> 3rd Qu.: 3.175 </td>
  </tr>
  <tr>
   <td style="text-align:left;">  </td>
   <td style="text-align:left;"> NA </td>
   <td style="text-align:left;"> Max.   :106 </td>
   <td style="text-align:left;"> Max.   :1223 </td>
   <td style="text-align:left;"> Max.   :126 </td>
   <td style="text-align:left;"> Max.   :1223 </td>
   <td style="text-align:left;"> Max.   :90953 </td>
   <td style="text-align:left;"> Max.   :442.44 </td>
   <td style="text-align:left;"> Max.   :0.37088 </td>
   <td style="text-align:left;"> Max.   :1.9596 </td>
   <td style="text-align:left;"> Max.   :0 </td>
   <td style="text-align:left;"> Max.   :43.44 </td>
   <td style="text-align:left;"> Max.   :4.0352 </td>
   <td style="text-align:left;"> Max.   :114577 </td>
   <td style="text-align:left;"> Max.   :114577 </td>
   <td style="text-align:left;"> Max.   :2427.94 </td>
   <td style="text-align:left;"> Max.   :441.14 </td>
   <td style="text-align:left;"> Max.   :8.597 </td>
   <td style="text-align:left;"> Max.   :7.483 </td>
   <td style="text-align:left;"> Max.   :43.44 </td>
   <td style="text-align:left;"> Max.   :173.3 </td>
   <td style="text-align:left;"> Max.   :9.5676 </td>
   <td style="text-align:left;"> Max.   :28.0 </td>
   <td style="text-align:left;"> Max.   :2.3e+09 </td>
   <td style="text-align:left;"> Max.   :4.0e+17 </td>
   <td style="text-align:left;"> Max.   :2.9e+25 </td>
   <td style="text-align:left;"> Max.   :4.6e+21 </td>
   <td style="text-align:left;"> Max.   :33.57 </td>
   <td style="text-align:left;"> Max.   :6.011 </td>
   <td style="text-align:left;"> Max.   :0.41105 </td>
   <td style="text-align:left;"> Max.   :1.6e+02 </td>
   <td style="text-align:left;"> Max.   :1.6e+14 </td>
   <td style="text-align:left;"> Max.   :2.6e+25 </td>
   <td style="text-align:left;"> Max.   :2.7e+28 </td>
   <td style="text-align:left;"> Max.   :2.8e+21 </td>
   <td style="text-align:left;"> Max.   :1.5e+21 </td>
   <td style="text-align:left;"> Max.   :3.7e+23 </td>
   <td style="text-align:left;"> Max.   :2.0e+03 </td>
   <td style="text-align:left;"> Max.   :1.9e+02 </td>
   <td style="text-align:left;"> Max.   :3892018 </td>
   <td style="text-align:left;"> Max.   :1.6e+02 </td>
   <td style="text-align:left;"> Max.   :1.8e+02 </td>
   <td style="text-align:left;"> Max.   :6.6e+04 </td>
   <td style="text-align:left;"> Max.   :303493 </td>
   <td style="text-align:left;"> Max.   : 7.0e+01 </td>
   <td style="text-align:left;"> Max.   :0.99422 </td>
   <td style="text-align:left;"> Max.   :1.00000 </td>
   <td style="text-align:left;"> Max.   :1.000 </td>
   <td style="text-align:left;"> Max.   :202.8 </td>
   <td style="text-align:left;"> Max.   :34.5 </td>
   <td style="text-align:left;"> Max.   :19.93 </td>
   <td style="text-align:left;"> Max.   :3.7e+08 </td>
   <td style="text-align:left;"> Max.   :8.7e+15 </td>
   <td style="text-align:left;"> Max.   :4.2e+22 </td>
   <td style="text-align:left;"> Max.   :2.0e+20 </td>
   <td style="text-align:left;"> Max.   :186.777 </td>
   <td style="text-align:left;"> Max.   :9.53689 </td>
   <td style="text-align:left;"> Max.   :1.8e+00 </td>
   <td style="text-align:left;"> Max.   :1.3e+03 </td>
   <td style="text-align:left;"> Max.   :1.3e+13 </td>
   <td style="text-align:left;"> Max.   :7.8e+23 </td>
   <td style="text-align:left;"> Max.   :1.6e+26 </td>
   <td style="text-align:left;"> Max.   :1.1e+20 </td>
   <td style="text-align:left;"> Max.   :4.5e+19 </td>
   <td style="text-align:left;"> Max.   :2.3e+21 </td>
   <td style="text-align:left;"> Max.   :6.1e+04 </td>
   <td style="text-align:left;"> Max.   :2.2e+03 </td>
   <td style="text-align:left;"> Max.   :3.7e+09 </td>
   <td style="text-align:left;"> Max.   :1.3e+03 </td>
   <td style="text-align:left;"> Max.   :1.3e+03 </td>
   <td style="text-align:left;"> Max.   :1.1e+07 </td>
   <td style="text-align:left;"> Max.   :55142 </td>
   <td style="text-align:left;"> Max.   : 9.0e+01 </td>
   <td style="text-align:left;"> Max.   :1.0e+00 </td>
   <td style="text-align:left;"> Max.   :1.00000 </td>
   <td style="text-align:left;"> Max.   :1.000 </td>
   <td style="text-align:left;"> Max.   :202.5 </td>
   <td style="text-align:left;"> Max.   :32.8 </td>
   <td style="text-align:left;"> Max.   :19.38 </td>
   <td style="text-align:left;"> Max.   :558.7 </td>
   <td style="text-align:left;"> Max.   :269.2 </td>
   <td style="text-align:left;"> Max.   :366.99 </td>
   <td style="text-align:left;"> Max.   :446.1 </td>
   <td style="text-align:left;"> Max.   :208.13 </td>
   <td style="text-align:left;"> Max.   :455.1 </td>
   <td style="text-align:left;"> Max.   :476.2 </td>
   <td style="text-align:left;"> Max.   :297.3 </td>
   <td style="text-align:left;"> Max.   :299.013 </td>
   <td style="text-align:left;"> Max.   :420.8 </td>
   <td style="text-align:left;"> Max.   :608.4 </td>
   <td style="text-align:left;"> Max.   :326.5 </td>
   <td style="text-align:left;"> Max.   :562 </td>
   <td style="text-align:left;"> Max.   :315.57 </td>
   <td style="text-align:left;"> Max.   :407.5 </td>
   <td style="text-align:left;"> Max.   :534.5 </td>
   <td style="text-align:left;"> Max.   :2.24 </td>
   <td style="text-align:left;"> Max.   :465.6 </td>
   <td style="text-align:left;"> Max.   :530.3 </td>
   <td style="text-align:left;"> Max.   :313.09 </td>
   <td style="text-align:left;"> Max.   :204.7 </td>
   <td style="text-align:left;"> Max.   :114.7 </td>
   <td style="text-align:left;"> Max.   :121.89 </td>
   <td style="text-align:left;"> Max.   :155.7 </td>
   <td style="text-align:left;"> Max.   :88.54 </td>
   <td style="text-align:left;"> Max.   :171.2 </td>
   <td style="text-align:left;"> Max.   :165.8 </td>
   <td style="text-align:left;"> Max.   :120.9 </td>
   <td style="text-align:left;"> Max.   :1.2e+02 </td>
   <td style="text-align:left;"> Max.   :176.3 </td>
   <td style="text-align:left;"> Max.   :260.1 </td>
   <td style="text-align:left;"> Max.   :135.3 </td>
   <td style="text-align:left;"> Max.   :228.2 </td>
   <td style="text-align:left;"> Max.   :118.23 </td>
   <td style="text-align:left;"> Max.   :156.8 </td>
   <td style="text-align:left;"> Max.   :228 </td>
   <td style="text-align:left;"> Max.   :2.25 </td>
   <td style="text-align:left;"> Max.   :181.06 </td>
   <td style="text-align:left;"> Max.   :186.7 </td>
   <td style="text-align:left;"> Max.   :1.2e+02 </td>
   <td style="text-align:left;"> Max.   :69202 </td>
   <td style="text-align:left;"> Max.   :69202 </td>
   <td style="text-align:left;"> Max.   :1996.25 </td>
   <td style="text-align:left;"> Max.   :395.695 </td>
   <td style="text-align:left;"> Max.   :8.857 </td>
   <td style="text-align:left;"> Max.   :7.7455 </td>
   <td style="text-align:left;"> Max.   :43.44 </td>
   <td style="text-align:left;"> Max.   :173.3 </td>
   <td style="text-align:left;"> Max.   :9.7427 </td>
   <td style="text-align:left;"> Max.   :24.0 </td>
   <td style="text-align:left;"> Max.   :1.8e+09 </td>
   <td style="text-align:left;"> Max.   :1.4e+17 </td>
   <td style="text-align:left;"> Max.   :5.9e+24 </td>
   <td style="text-align:left;"> Max.   :1.9e+21 </td>
   <td style="text-align:left;"> Max.   :38.43 </td>
   <td style="text-align:left;"> Max.   :11.751 </td>
   <td style="text-align:left;"> Max.   :1.23161 </td>
   <td style="text-align:left;"> Max.   :5.5e+02 </td>
   <td style="text-align:left;"> Max.   :1.3e+14 </td>
   <td style="text-align:left;"> Max.   :8.0e+24 </td>
   <td style="text-align:left;"> Max.   :1.7e+28 </td>
   <td style="text-align:left;"> Max.   :9.6e+20 </td>
   <td style="text-align:left;"> Max.   :4.4e+20 </td>
   <td style="text-align:left;"> Max.   :2.4e+23 </td>
   <td style="text-align:left;"> Max.   :2.6e+03 </td>
   <td style="text-align:left;"> Max.   :1.0e+03 </td>
   <td style="text-align:left;"> Max.   :6717223 </td>
   <td style="text-align:left;"> Max.   :6.1e+02 </td>
   <td style="text-align:left;"> Max.   :640.948 </td>
   <td style="text-align:left;"> Max.   :1.0e+05 </td>
   <td style="text-align:left;"> Max.   :249531 </td>
   <td style="text-align:left;"> Max.   : 6.8e+01 </td>
   <td style="text-align:left;"> Max.   :9.9e-01 </td>
   <td style="text-align:left;"> Max.   :1.00000 </td>
   <td style="text-align:left;"> Max.   :1.000 </td>
   <td style="text-align:left;"> Max.   :202.4 </td>
   <td style="text-align:left;"> Max.   :32.06 </td>
   <td style="text-align:left;"> Max.   :19.06 </td>
   <td style="text-align:left;"> Max.   :3.3e+08 </td>
   <td style="text-align:left;"> Max.   :6.8e+15 </td>
   <td style="text-align:left;"> Max.   :2.7e+22 </td>
   <td style="text-align:left;"> Max.   :1.5e+20 </td>
   <td style="text-align:left;"> Max.   :205.735 </td>
   <td style="text-align:left;"> Max.   :1.4e+01 </td>
   <td style="text-align:left;"> Max.   :3.1e+00 </td>
   <td style="text-align:left;"> Max.   :4.7e+03 </td>
   <td style="text-align:left;"> Max.   :1.1e+13 </td>
   <td style="text-align:left;"> Max.   :6.2e+23 </td>
   <td style="text-align:left;"> Max.   :1.2e+26 </td>
   <td style="text-align:left;"> Max.   :8.4e+19 </td>
   <td style="text-align:left;"> Max.   :3.8e+19 </td>
   <td style="text-align:left;"> Max.   :1.7e+21 </td>
   <td style="text-align:left;"> Max.   :7.4e+04 </td>
   <td style="text-align:left;"> Max.   :3.0e+04 </td>
   <td style="text-align:left;"> Max.   :5.5e+09 </td>
   <td style="text-align:left;"> Max.   :4.7e+03 </td>
   <td style="text-align:left;"> Max.   :4.7e+03 </td>
   <td style="text-align:left;"> Max.   :1.5e+07 </td>
   <td style="text-align:left;"> Max.   :49462 </td>
   <td style="text-align:left;"> Max.   : 1.0e+02 </td>
   <td style="text-align:left;"> Max.   :1.0e+00 </td>
   <td style="text-align:left;"> Max.   :1.00000 </td>
   <td style="text-align:left;"> Max.   :1.000 </td>
   <td style="text-align:left;"> Max.   :202.2 </td>
   <td style="text-align:left;"> Max.   :30.77 </td>
   <td style="text-align:left;"> Max.   :18.69 </td>
   <td style="text-align:left;"> Max.   :470.3 </td>
   <td style="text-align:left;"> Max.   :244 </td>
   <td style="text-align:left;"> Max.   :291.56 </td>
   <td style="text-align:left;"> Max.   :371.0 </td>
   <td style="text-align:left;"> Max.   :191.76 </td>
   <td style="text-align:left;"> Max.   :403.3 </td>
   <td style="text-align:left;"> Max.   :371.7 </td>
   <td style="text-align:left;"> Max.   :260.7 </td>
   <td style="text-align:left;"> Max.   :263.618 </td>
   <td style="text-align:left;"> Max.   :346.4 </td>
   <td style="text-align:left;"> Max.   :483.6 </td>
   <td style="text-align:left;"> Max.   :284.60 </td>
   <td style="text-align:left;"> Max.   :451.6 </td>
   <td style="text-align:left;"> Max.   :295.23 </td>
   <td style="text-align:left;"> Max.   :368.1 </td>
   <td style="text-align:left;"> Max.   :422.5 </td>
   <td style="text-align:left;"> Max.   :3.2 </td>
   <td style="text-align:left;"> Max.   :375 </td>
   <td style="text-align:left;"> Max.   :447 </td>
   <td style="text-align:left;"> Max.   :275.172 </td>
   <td style="text-align:left;"> Max.   :195.0 </td>
   <td style="text-align:left;"> Max.   :108.7 </td>
   <td style="text-align:left;"> Max.   :111.8 </td>
   <td style="text-align:left;"> Max.   :147.5 </td>
   <td style="text-align:left;"> Max.   :86.32 </td>
   <td style="text-align:left;"> Max.   :163.5 </td>
   <td style="text-align:left;"> Max.   :158.59 </td>
   <td style="text-align:left;"> Max.   :117.4 </td>
   <td style="text-align:left;"> Max.   :117.454 </td>
   <td style="text-align:left;"> Max.   :170.9 </td>
   <td style="text-align:left;"> Max.   :252.0 </td>
   <td style="text-align:left;"> Max.   :128.13 </td>
   <td style="text-align:left;"> Max.   :218.4 </td>
   <td style="text-align:left;"> Max.   :112.8 </td>
   <td style="text-align:left;"> Max.   :147.4 </td>
   <td style="text-align:left;"> Max.   :221.7 </td>
   <td style="text-align:left;"> Max.   :3.25 </td>
   <td style="text-align:left;"> Max.   :175 </td>
   <td style="text-align:left;"> Max.   :182.2 </td>
   <td style="text-align:left;"> Max.   :1.2e+02 </td>
   <td style="text-align:left;"> Max.   :45564 </td>
   <td style="text-align:left;"> Max.   :45564 </td>
   <td style="text-align:left;"> Max.   :1632.54 </td>
   <td style="text-align:left;"> Max.   :351.187 </td>
   <td style="text-align:left;"> Max.   :8.973 </td>
   <td style="text-align:left;"> Max.   :7.95879 </td>
   <td style="text-align:left;"> Max.   :43.44 </td>
   <td style="text-align:left;"> Max.   :173.3 </td>
   <td style="text-align:left;"> Max.   :9.94605 </td>
   <td style="text-align:left;"> Max.   :26.0 </td>
   <td style="text-align:left;"> Max.   :1.5e+09 </td>
   <td style="text-align:left;"> Max.   :5.9e+16 </td>
   <td style="text-align:left;"> Max.   :2.0e+24 </td>
   <td style="text-align:left;"> Max.   : 1.2e+21 </td>
   <td style="text-align:left;"> Max.   :68.88 </td>
   <td style="text-align:left;"> Max.   :21.693 </td>
   <td style="text-align:left;"> Max.   :6.41781 </td>
   <td style="text-align:left;"> Max.   :3.4e+03 </td>
   <td style="text-align:left;"> Max.   :1.1e+14 </td>
   <td style="text-align:left;"> Max.   :5.3e+24 </td>
   <td style="text-align:left;"> Max.   :1.1e+28 </td>
   <td style="text-align:left;"> Max.   : 6.5e+20 </td>
   <td style="text-align:left;"> Max.   :2.6e+20 </td>
   <td style="text-align:left;"> Max.   :1.6e+23 </td>
   <td style="text-align:left;"> Max.   :8.4e+03 </td>
   <td style="text-align:left;"> Max.   :4.8e+04 </td>
   <td style="text-align:left;"> Max.   :7.0e+07 </td>
   <td style="text-align:left;"> Max.   :3.4e+03 </td>
   <td style="text-align:left;"> Max.   :3.4e+03 </td>
   <td style="text-align:left;"> Max.   :5.8e+05 </td>
   <td style="text-align:left;"> Max.   :204067 </td>
   <td style="text-align:left;"> Max.   : 1.1e+02 </td>
   <td style="text-align:left;"> Max.   :1.0e+00 </td>
   <td style="text-align:left;"> Max.   :1.0e+00 </td>
   <td style="text-align:left;"> Max.   :1.00 </td>
   <td style="text-align:left;"> Max.   :201.81 </td>
   <td style="text-align:left;"> Max.   :29.65 </td>
   <td style="text-align:left;"> Max.   :18.22 </td>
   <td style="text-align:left;"> Max.   :3.0e+08 </td>
   <td style="text-align:left;"> Max.   :5.2e+15 </td>
   <td style="text-align:left;"> Max.   :1.7e+22 </td>
   <td style="text-align:left;"> Max.   : 1.2e+20 </td>
   <td style="text-align:left;"> Max.   :382.977 </td>
   <td style="text-align:left;"> Max.   :1.2e+02 </td>
   <td style="text-align:left;"> Max.   :1.1e+02 </td>
   <td style="text-align:left;"> Max.   :2.0e+05 </td>
   <td style="text-align:left;"> Max.   :9.4e+12 </td>
   <td style="text-align:left;"> Max.   :4.8e+23 </td>
   <td style="text-align:left;"> Max.   :8.8e+25 </td>
   <td style="text-align:left;"> Max.   : 6.7e+19 </td>
   <td style="text-align:left;"> Max.   :3.2e+19 </td>
   <td style="text-align:left;"> Max.   :1.2e+21 </td>
   <td style="text-align:left;"> Max.   :2.6e+05 </td>
   <td style="text-align:left;"> Max.   :2.2e+06 </td>
   <td style="text-align:left;"> Max.   :6.7e+10 </td>
   <td style="text-align:left;"> Max.   :2.6e+05 </td>
   <td style="text-align:left;"> Max.   :3.1e+05 </td>
   <td style="text-align:left;"> Max.   :9.9e+07 </td>
   <td style="text-align:left;"> Max.   :43898.3 </td>
   <td style="text-align:left;"> Max.   : 1.2e+02 </td>
   <td style="text-align:left;"> Max.   :9.9e-01 </td>
   <td style="text-align:left;"> Max.   :1.0e+00 </td>
   <td style="text-align:left;"> Max.   :0.9999 </td>
   <td style="text-align:left;"> Max.   :201.71 </td>
   <td style="text-align:left;"> Max.   :28.75 </td>
   <td style="text-align:left;"> Max.   :17.68 </td>
   <td style="text-align:left;"> Max.   :414.4 </td>
   <td style="text-align:left;"> Max.   :220.7 </td>
   <td style="text-align:left;"> Max.   :228.84 </td>
   <td style="text-align:left;"> Max.   :319.7 </td>
   <td style="text-align:left;"> Max.   :195.16 </td>
   <td style="text-align:left;"> Max.   :345 </td>
   <td style="text-align:left;"> Max.   :322.34 </td>
   <td style="text-align:left;"> Max.   :262.0 </td>
   <td style="text-align:left;"> Max.   :1.9e+02 </td>
   <td style="text-align:left;"> Max.   :315.9 </td>
   <td style="text-align:left;"> Max.   :441.0 </td>
   <td style="text-align:left;"> Max.   :256.750 </td>
   <td style="text-align:left;"> Max.   :402.7 </td>
   <td style="text-align:left;"> Max.   :241.83 </td>
   <td style="text-align:left;"> Max.   :300.4 </td>
   <td style="text-align:left;"> Max.   :379.8 </td>
   <td style="text-align:left;"> Max.   :4.57 </td>
   <td style="text-align:left;"> Max.   :345.04 </td>
   <td style="text-align:left;"> Max.   :377.5 </td>
   <td style="text-align:left;"> Max.   :2.1e+02 </td>
   <td style="text-align:left;"> Max.   :189.1 </td>
   <td style="text-align:left;"> Max.   :102.37 </td>
   <td style="text-align:left;"> Max.   :107.1 </td>
   <td style="text-align:left;"> Max.   :141.5 </td>
   <td style="text-align:left;"> Max.   :83.38 </td>
   <td style="text-align:left;"> Max.   :156.1 </td>
   <td style="text-align:left;"> Max.   :153.30 </td>
   <td style="text-align:left;"> Max.   :113.0 </td>
   <td style="text-align:left;"> Max.   :1.2e+02 </td>
   <td style="text-align:left;"> Max.   :165.28 </td>
   <td style="text-align:left;"> Max.   :243.0 </td>
   <td style="text-align:left;"> Max.   :120.557 </td>
   <td style="text-align:left;"> Max.   :208.21 </td>
   <td style="text-align:left;"> Max.   :105.93 </td>
   <td style="text-align:left;"> Max.   :136.5 </td>
   <td style="text-align:left;"> Max.   :214.21 </td>
   <td style="text-align:left;"> Max.   :4.58 </td>
   <td style="text-align:left;"> Max.   :168.9 </td>
   <td style="text-align:left;"> Max.   :176.8 </td>
   <td style="text-align:left;"> Max.   :1.2e+02 </td>
   <td style="text-align:left;"> Max.   :7.60 </td>
   <td style="text-align:left;"> Max.   : 2.6e-07 </td>
   <td style="text-align:left;"> Max.   :0.797 </td>
   <td style="text-align:left;"> Max.   :0.6360 </td>
   <td style="text-align:left;"> Max.   :-0.043 </td>
   <td style="text-align:left;"> Max.   :43.352 </td>
  </tr>
</tbody>
</table></div>

#Ograniczenie zakresu

Poni¿szy kod ogranicza liczbê wierszy, pozostawiaj¹c tylko 50 najczêœciej wystêpuj¹cych wartoœci atrybutu res_name.


```r
top_50 <- data_final %>% group_by(res_name) %>%
  summarize(ilosc=n()) %>% arrange(desc(ilosc)) %>% head(50) 
data_final_50 <- data_final %>% filter(res_name %in% top_50[['res_name']])
```

#Korelacja miêdzy zmiennymi

Poni¿ej przedstawiono korelacjê pomiêdzy niektórymi zmiennymi.


```r
cor_grouped <- data_final_50 %>% select(local_res_atom_non_h_count:local_skewness, -local_min)

cor_data <- melt(cor(cor_grouped))
ggplot(cor_data, aes(x = Var1, y = Var2, fill = value)) + geom_tile() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
```

![](ZED_files/figure-html/korelacja-1.png)<!-- -->

#Liczebnoœæ 50 najpopularniejszych klas


```r
knitr::kable(top_50, align = 'l') %>%
  kable_styling("striped", full_width = FALSE) %>% 
  column_spec(1, width = "2.5cm")
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> res_name </th>
   <th style="text-align:left;"> ilosc </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> SO4 </td>
   <td style="text-align:left;"> 28354 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> GOL </td>
   <td style="text-align:left;"> 19163 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> EDO </td>
   <td style="text-align:left;"> 14530 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> CL </td>
   <td style="text-align:left;"> 11443 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAG </td>
   <td style="text-align:left;"> 11432 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> CA </td>
   <td style="text-align:left;"> 10528 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> ZN </td>
   <td style="text-align:left;"> 9989 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MG </td>
   <td style="text-align:left;"> 6957 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> HEM </td>
   <td style="text-align:left;"> 5305 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> PO4 </td>
   <td style="text-align:left;"> 5287 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> ACT </td>
   <td style="text-align:left;"> 3769 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> IOD </td>
   <td style="text-align:left;"> 3374 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> DMS </td>
   <td style="text-align:left;"> 3364 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAD </td>
   <td style="text-align:left;"> 2556 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> K </td>
   <td style="text-align:left;"> 2317 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> FAD </td>
   <td style="text-align:left;"> 2248 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> PEG </td>
   <td style="text-align:left;"> 2221 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MN </td>
   <td style="text-align:left;"> 2075 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> ADP </td>
   <td style="text-align:left;"> 1822 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAP </td>
   <td style="text-align:left;"> 1720 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> CD </td>
   <td style="text-align:left;"> 1696 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MLY </td>
   <td style="text-align:left;"> 1676 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> UNX </td>
   <td style="text-align:left;"> 1479 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MPD </td>
   <td style="text-align:left;"> 1380 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> FMT </td>
   <td style="text-align:left;"> 1328 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MES </td>
   <td style="text-align:left;"> 1326 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> PG4 </td>
   <td style="text-align:left;"> 1325 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> CU </td>
   <td style="text-align:left;"> 1140 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MAN </td>
   <td style="text-align:left;"> 1106 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> BR </td>
   <td style="text-align:left;"> 1101 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> ATP </td>
   <td style="text-align:left;"> 1056 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> FMN </td>
   <td style="text-align:left;"> 1050 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> COA </td>
   <td style="text-align:left;"> 1049 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> 1PE </td>
   <td style="text-align:left;"> 1020 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> EPE </td>
   <td style="text-align:left;"> 965 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> CLA </td>
   <td style="text-align:left;"> 939 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NDP </td>
   <td style="text-align:left;"> 912 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NI </td>
   <td style="text-align:left;"> 898 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NO3 </td>
   <td style="text-align:left;"> 880 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> SF4 </td>
   <td style="text-align:left;"> 863 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> ACY </td>
   <td style="text-align:left;"> 813 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> TRS </td>
   <td style="text-align:left;"> 805 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> SAH </td>
   <td style="text-align:left;"> 801 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> GDP </td>
   <td style="text-align:left;"> 798 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> PGE </td>
   <td style="text-align:left;"> 760 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> PLP </td>
   <td style="text-align:left;"> 751 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> HEC </td>
   <td style="text-align:left;"> 747 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> FE </td>
   <td style="text-align:left;"> 734 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> CIT </td>
   <td style="text-align:left;"> 719 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> FE2 </td>
   <td style="text-align:left;"> 694 </td>
  </tr>
</tbody>
</table>

#Wykres rozk³adu liczby atomów


```r
plot_atom <- ggplot(data_final_50, aes(local_res_atom_non_h_count, fill = res_name)) + geom_histogram()
ggplotly(plot_atom)
```

<!--html_preserve--><div id="htmlwidget-49c968e2f00b87e48d74" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-49c968e2f00b87e48d74">{"x":{"data":[{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[54433,1375,57080,21503,3309,4725,12945,1902,21,4,26,28,2719,798,2191,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,4,12,65,71,231,76,561,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:    12<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:    65<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:    71<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:   231<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:    76<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:   561<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: 1PE"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"1PE","legendgroup":"1PE","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[54433,1367,53319,21503,3309,4725,12945,1902,21,4,26,28,2719,798,2191,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,8,3761,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     8<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:  3761<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACT"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(243,123,89,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ACT","legendgroup":"ACT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[54433,1358,52515,21503,3309,4725,12945,1902,21,4,26,28,2719,798,2191,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,9,804,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     9<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:   804<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ACY"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(237,129,65,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ACY","legendgroup":"ACY","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[54433,1358,52515,21503,3305,4725,12945,1902,16,4,26,26,908,798,2191,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,4,0,0,0,5,0,0,2,1811,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     5<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:  1811<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ADP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(231,134,27,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ADP","legendgroup":"ADP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[54433,1358,52515,21503,3303,4723,12933,1902,16,3,24,26,901,795,1164,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,2,2,12,0,0,1,2,0,7,3,1027,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:    12<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     7<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     3<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:  1027<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: ATP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(224,139,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ATP","legendgroup":"ATP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[53332,1358,52515,21503,3303,4723,12933,1902,16,3,24,26,901,795,1164,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[1101,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  1101<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: BR"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(216,144,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"BR","legendgroup":"BR","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[42804,1358,52515,21503,3303,4723,12933,1902,16,3,24,26,901,795,1164,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[10528,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count: 10528<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CA"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(207,148,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CA","legendgroup":"CA","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[41108,1358,52515,21503,3303,4723,12933,1902,16,3,24,26,901,795,1164,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[1696,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  1696<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(197,153,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CD","legendgroup":"CD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[41108,1358,52515,21503,3303,4719,12218,1902,16,3,24,26,901,795,1164,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,4,715,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:   715<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CIT"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(187,157,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CIT","legendgroup":"CIT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[29665,1358,52515,21503,3303,4719,12218,1902,16,3,24,26,901,795,1164,3,62,17,70,6063,2457,28,3395,32,2241,41,5,13,4,755],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[11443,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count: 11443<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CL"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(175,161,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CL","legendgroup":"CL","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[29665,1358,52515,21503,3303,4719,12218,1902,16,3,24,6,901,795,1164,3,55,17,70,6059,2448,1,3391,0,2223,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,0,0,0,0,0,20,0,0,0,0,7,0,0,4,9,27,4,32,18,41,5,13,4,755],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:    20<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     7<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     9<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:    27<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:    32<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:    18<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:    41<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     5<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:    13<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: CLA","count:   755<br />local_res_atom_non_h_count: NA<br />res_name: CLA"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(163,165,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CLA","legendgroup":"CLA","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[29665,1358,52515,21503,3302,4719,12218,1902,4,2,14,6,892,795,1137,1,55,10,65,6042,2447,1,2434,0,2223,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,1,0,0,0,12,1,10,0,9,0,27,2,0,7,5,17,1,0,957,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:    12<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:    10<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     9<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:    27<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     7<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     5<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:    17<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:   957<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: COA"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(149,169,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"COA","legendgroup":"COA","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[28525,1358,52515,21503,3302,4719,12218,1902,4,2,14,6,892,795,1137,1,55,10,65,6042,2447,1,2434,0,2223,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[1140,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  1140<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: CU"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(133,173,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CU","legendgroup":"CU","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[28525,1356,49153,21503,3302,4719,12218,1902,4,2,14,6,892,795,1137,1,55,10,65,6042,2447,1,2434,0,2223,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,2,3362,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:  3362<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: DMS"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(114,176,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"DMS","legendgroup":"DMS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[28525,1350,34629,21503,3302,4719,12218,1902,4,2,14,6,892,795,1137,1,55,10,65,6042,2447,1,2434,0,2223,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,6,14524,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     6<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count: 14524<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EDO"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(91,179,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"EDO","legendgroup":"EDO","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[28525,1350,34608,21495,3295,4684,12204,1022,4,2,14,6,892,795,1137,1,55,10,65,6042,2447,1,2434,0,2223,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,21,8,7,35,14,880,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:    21<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     8<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     7<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:    35<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:    14<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:   880<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: EPE"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(57,182,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"EPE","legendgroup":"EPE","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[28525,1350,34608,21495,3295,4684,12204,1022,4,2,14,6,892,795,1137,1,34,8,65,6040,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,21,2,0,2,0,0,0,0,2223,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:    21<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:  2223<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FAD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,184,31,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FAD","legendgroup":"FAD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[27791,1350,34608,21495,3295,4684,12204,1022,4,2,14,6,892,795,1137,1,34,8,65,6040,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[734,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:   734<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,186,66,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FE","legendgroup":"FE","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[27097,1350,34608,21495,3295,4684,12204,1022,4,2,14,6,892,795,1137,1,34,8,65,6040,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[694,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:   694<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: FE2"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,188,89,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FE2","legendgroup":"FE2","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[27097,1350,34608,21495,3295,4684,12204,1022,4,1,14,4,892,795,90,1,34,8,65,6040,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,0,0,0,1,0,2,0,0,1047,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:  1047<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMN"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,190,108,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FMN","legendgroup":"FMN","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[27097,22,34608,21495,3295,4684,12204,1022,4,1,14,4,892,795,90,1,34,8,65,6040,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,1328,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:  1328<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: FMT"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,191,125,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FMT","legendgroup":"FMT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[27097,22,34608,21495,3295,4683,12204,1022,2,1,14,3,892,1,90,1,34,8,65,6040,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,1,0,0,2,0,0,1,0,794,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:   794<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GDP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,192,141,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"GDP","legendgroup":"GDP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[27097,21,34581,2360,3295,4683,12204,1022,2,1,14,3,892,1,90,1,34,8,65,6040,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,1,27,19135,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:    27<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count: 19135<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: GOL"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,193,156,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"GOL","legendgroup":"GOL","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[27097,21,34581,2360,3295,4683,12204,1022,2,1,14,3,892,1,90,1,34,0,60,5306,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,5,734,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     8<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     5<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:   734<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEC"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,193,170,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"HEC","legendgroup":"HEC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[27097,21,34581,2360,3295,4683,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,5303,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:  5303<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: HEM"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,192,184,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"HEM","legendgroup":"HEM","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[23723,21,34581,2360,3295,4683,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[3374,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  3374<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: IOD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,191,196,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"IOD","legendgroup":"IOD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[21406,21,34581,2360,3295,4683,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[2317,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  2317<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: K"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,189,208,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"K","legendgroup":"K","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[21406,21,34581,2360,3293,3579,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,2,1104,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:  1104<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MAN"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,187,219,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MAN","legendgroup":"MAN","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[21406,21,34577,2354,3291,2265,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,4,6,2,1314,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     6<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:  1314<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MES"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,184,229,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MES","legendgroup":"MES","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[14449,21,34577,2354,3291,2265,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[6957,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  6957<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MG"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,180,239,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MG","legendgroup":"MG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[14449,21,34556,2333,3069,853,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,21,21,222,1412,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:    21<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:    21<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:   222<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:  1412<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MLY"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,176,246,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MLY","legendgroup":"MLY","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[12374,21,34556,2333,3069,853,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[2075,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  2075<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: MN"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,171,253,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MN","legendgroup":"MN","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[12374,21,34554,2329,1695,853,12204,1022,2,1,14,3,892,1,90,0,34,0,59,3,2447,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,2,4,1374,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:  1374<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: MPD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,165,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MPD","legendgroup":"MPD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[12374,21,34554,2329,1695,846,12204,1022,2,1,9,3,823,0,90,0,0,0,59,0,10,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,7,0,0,0,0,5,0,69,1,0,0,34,0,0,3,2437,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     7<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     5<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:    69<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:    34<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     3<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:  2437<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(82,158,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NAD","legendgroup":"NAD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[12373,21,34554,2329,1694,843,1059,740,2,1,9,3,823,0,90,0,0,0,59,0,10,1,2434,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[1,0,0,0,1,3,11145,282,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     1<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     3<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count: 11145<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:   282<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAG"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(121,151,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NAG","legendgroup":"NAG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[12373,21,34554,2329,1694,843,1059,740,0,1,3,3,807,0,21,0,0,0,10,0,0,1,866,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,0,0,2,0,6,0,16,0,69,0,0,0,49,0,10,0,1568,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     6<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:    16<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:    69<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:    49<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:    10<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:  1568<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NAP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(149,144,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NAP","legendgroup":"NAP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[12373,21,34554,2329,1694,843,1059,740,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,0,0,0,0,3,3,8,0,21,0,0,0,10,0,0,1,866,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     3<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     3<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     8<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:    21<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:    10<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:   866<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: NDP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(172,136,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NDP","legendgroup":"NDP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11475,21,34554,2329,1694,843,1059,740,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[898,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:   898<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: NI"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(191,128,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NI","legendgroup":"NI","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11475,21,33674,2329,1694,843,1059,740,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,880,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:   880<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3","count:     0<br />local_res_atom_non_h_count:  4<br />res_name: NO3"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(207,120,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NO3","legendgroup":"NO3","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11475,19,33651,133,1694,843,1059,740,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,2,23,2196,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     2<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:    23<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:  2196<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PEG"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(220,113,250,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PEG","legendgroup":"PEG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11475,19,33637,31,1660,715,12,740,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,14,102,34,128,1047,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:    14<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:   102<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:    34<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:   128<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:  1047<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PG4"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(231,107,243,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PG4","legendgroup":"PG4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11475,19,33620,9,1654,0,12,740,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,17,22,6,715,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:    17<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:    22<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     6<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:   715<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PGE"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(240,102,234,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PGE","legendgroup":"PGE","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11475,19,33620,9,1654,0,0,1,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,12,739,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:    12<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:   739<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PLP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(247,99,224,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PLP","legendgroup":"PLP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11473,12,28342,9,1654,0,0,1,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[2,7,5278,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     2<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     7<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:  5278<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: PO4"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(252,97,213,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PO4","legendgroup":"PO4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11473,12,28342,9,1654,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,0,0,0,0,0,1,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:   799<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SAH"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,97,201,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"SAH","legendgroup":"SAH","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11473,12,28341,6,795,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,1,3,859,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     1<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     3<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:   859<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SF4"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,98,188,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"SF4","legendgroup":"SF4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11468,0,4,6,795,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[5,12,28337,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     5<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:    12<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count: 28337<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: SO4"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,101,174,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"SO4","legendgroup":"SO4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[11468,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[0,0,4,6,795,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     4<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     6<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:   795<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_count: NA<br />res_name: TRS"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,104,159,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"TRS","legendgroup":"TRS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[9989,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[1479,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  1479<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: UNX"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,108,144,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"UNX","legendgroup":"UNX","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172413,2.20689655172413,2.20689655172414,2.20689655172416],"base":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,2.20689655172414,4.41379310344828,6.62068965517241,8.82758620689655,11.0344827586207,13.2413793103448,15.448275862069,17.6551724137931,19.8620689655172,22.0689655172414,24.2758620689655,26.4827586206897,28.6896551724138,30.8965517241379,33.1034482758621,35.3103448275862,37.5172413793103,39.7241379310345,41.9310344827586,44.1379310344828,46.3448275862069,48.551724137931,50.7586206896552,52.9655172413793,55.1724137931034,57.3793103448276,59.5862068965517,61.7931034482759,64],"y":[9989,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  9989<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN","count:     0<br />local_res_atom_non_h_count:  1<br />res_name: ZN"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(252,113,127,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ZN","legendgroup":"ZN","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":26.2283105022831,"r":7.30593607305936,"b":40.1826484018265,"l":54.7945205479452},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-4.41379310344828,68.4137931034483],"tickmode":"array","ticktext":["0","20","40","60"],"tickvals":[0,20,40,60],"categoryorder":"array","categoryarray":["0","20","40","60"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"y","title":"local_res_atom_non_h_count","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-2854.6,59946.6],"tickmode":"array","ticktext":["0","20000","40000"],"tickvals":[0,20000,40000],"categoryorder":"array","categoryarray":["0","20000","40000"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"x","title":"count","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"y":0.93503937007874},"annotations":[{"text":"res_name","x":1.02,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"left","yanchor":"bottom","legendTitle":true}],"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":[{"name":"Collaborate","icon":{"width":1000,"ascent":500,"descent":-50,"path":"M487 375c7-10 9-23 5-36l-79-259c-3-12-11-23-22-31-11-8-22-12-35-12l-263 0c-15 0-29 5-43 15-13 10-23 23-28 37-5 13-5 25-1 37 0 0 0 3 1 7 1 5 1 8 1 11 0 2 0 4-1 6 0 3-1 5-1 6 1 2 2 4 3 6 1 2 2 4 4 6 2 3 4 5 5 7 5 7 9 16 13 26 4 10 7 19 9 26 0 2 0 5 0 9-1 4-1 6 0 8 0 2 2 5 4 8 3 3 5 5 5 7 4 6 8 15 12 26 4 11 7 19 7 26 1 1 0 4 0 9-1 4-1 7 0 8 1 2 3 5 6 8 4 4 6 6 6 7 4 5 8 13 13 24 4 11 7 20 7 28 1 1 0 4 0 7-1 3-1 6-1 7 0 2 1 4 3 6 1 1 3 4 5 6 2 3 3 5 5 6 1 2 3 5 4 9 2 3 3 7 5 10 1 3 2 6 4 10 2 4 4 7 6 9 2 3 4 5 7 7 3 2 7 3 11 3 3 0 8 0 13-1l0-1c7 2 12 2 14 2l218 0c14 0 25-5 32-16 8-10 10-23 6-37l-79-259c-7-22-13-37-20-43-7-7-19-10-37-10l-248 0c-5 0-9-2-11-5-2-3-2-7 0-12 4-13 18-20 41-20l264 0c5 0 10 2 16 5 5 3 8 6 10 11l85 282c2 5 2 10 2 17 7-3 13-7 17-13z m-304 0c-1-3-1-5 0-7 1-1 3-2 6-2l174 0c2 0 4 1 7 2 2 2 4 4 5 7l6 18c0 3 0 5-1 7-1 1-3 2-6 2l-173 0c-3 0-5-1-8-2-2-2-4-4-4-7z m-24-73c-1-3-1-5 0-7 2-2 3-2 6-2l174 0c2 0 5 0 7 2 3 2 4 4 5 7l6 18c1 2 0 5-1 6-1 2-3 3-5 3l-174 0c-3 0-5-1-7-3-3-1-4-4-5-6z"},"click":"function(gd) { \n        // is this being viewed in RStudio?\n        if (location.search == '?viewer_pane=1') {\n          alert('To learn about plotly for collaboration, visit:\\n https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html');\n        } else {\n          window.open('https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html', '_blank');\n        }\n      }"}],"cloud":false},"source":"A","attrs":{"234822225c04":{"x":{},"fill":{},"type":"bar"}},"cur_data":"234822225c04","visdat":{"234822225c04":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"base_url":"https://plot.ly"},"evals":["config.modeBarButtonsToAdd.0.click"],"jsHooks":[]}</script><!--/html_preserve-->

#Wykres rozk³adu liczby elektronów


```r
plot_elect <- ggplot(data_final_50, aes(local_res_atom_non_h_electron_sum, fill = res_name)) + geom_histogram()
ggplotly(plot_elect)
```

<!--html_preserve--><div id="htmlwidget-f474ad921ae63e0a53f4" style="width:672px;height:480px;" class="plotly html-widget"></div>
<script type="application/json" data-for="htmlwidget-f474ad921ae63e0a53f4">{"x":{"data":[{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,31275,36960,61286,5825,3238,1252,13461,1624,26,9,27,872,810,10,2713,1055,11,1145,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,4,12,65,71,167,128,23,550,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     4<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:    12<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:    65<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:    71<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:   167<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:   128<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:    23<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:   550<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: 1PE"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(248,118,109,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"1PE","legendgroup":"1PE","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,31271,33195,61286,5825,3238,1252,13461,1624,26,9,27,872,810,10,2713,1055,11,1145,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,4,3765,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     4<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:  3765<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACT"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(243,123,89,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ACT","legendgroup":"ACT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,31262,32391,61286,5825,3238,1252,13461,1624,26,9,27,872,810,10,2713,1055,11,1145,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,9,804,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     9<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:   804<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ACY"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(237,129,65,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ACY","legendgroup":"ACY","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,31262,32391,61286,5825,3238,1248,13461,1621,26,8,26,872,810,8,902,1055,11,1145,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,4,0,3,0,1,1,0,0,2,1811,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     4<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:  1811<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ADP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(231,134,27,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ADP","legendgroup":"ADP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,31262,32391,61286,5825,3237,1246,13461,1620,13,8,26,870,810,8,895,1052,0,129,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,1,2,0,1,13,0,0,2,0,0,7,3,11,1016,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:    13<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     7<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:    11<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:  1016<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: ATP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(224,139,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ATP","legendgroup":"ATP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,31262,32391,60185,5825,3237,1246,13461,1620,13,8,26,870,810,8,895,1052,0,129,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,1101,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:  1101<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR","count:     0<br />local_res_atom_non_h_electron_sum: 35<br />res_name: BR"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(216,144,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"BR","legendgroup":"BR","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,20734,32391,60185,5825,3237,1246,13461,1620,13,8,26,870,810,8,895,1052,0,129,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,10528,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count: 10528<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA","count:     0<br />local_res_atom_non_h_electron_sum: 20<br />res_name: CA"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(207,148,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CA","legendgroup":"CA","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,20734,32391,58489,5825,3237,1246,13461,1620,13,8,26,870,810,8,895,1052,0,129,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,1696,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:  1696<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD","count:     0<br />local_res_atom_non_h_electron_sum: 48<br />res_name: CD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(197,153,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CD","legendgroup":"CD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,20734,32391,58489,5825,3235,1244,12746,1620,13,8,26,870,810,8,895,1052,0,129,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,2,2,715,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:   715<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CIT"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(187,157,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CIT","legendgroup":"CIT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9291,32391,58489,5825,3235,1244,12746,1620,13,8,26,870,810,8,895,1052,0,129,66,27,6086,43,2491,9,52,3404,2229,4,755],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,11443,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count: 11443<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL","count:     0<br />local_res_atom_non_h_electron_sum: 17<br />res_name: CL"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(175,161,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CL","legendgroup":"CL","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9291,32391,58489,5825,3235,1244,12746,1620,13,8,6,870,810,8,895,1045,0,129,62,25,6069,24,2457,3,0,3398,2216,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,0,0,0,0,20,0,0,0,0,7,0,0,4,2,17,19,34,6,52,6,13,4,755],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:    20<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     7<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     4<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:    17<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:    19<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:    34<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     6<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:    52<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     6<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:    13<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:     4<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA","count:   755<br />local_res_atom_non_h_electron_sum: NA<br />res_name: CLA"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(163,165,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CLA","legendgroup":"CLA","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9291,32391,58489,5825,3234,1244,12746,1620,1,8,3,870,802,8,886,1045,0,100,62,25,6057,18,2446,2,0,2441,2216,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,1,0,0,0,12,0,3,0,8,0,9,0,0,29,0,0,12,6,11,1,0,957,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:    12<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     8<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     9<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:    29<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:    12<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     6<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:    11<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:   957<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: COA"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(149,169,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"COA","legendgroup":"COA","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9291,31251,58489,5825,3234,1244,12746,1620,1,8,3,870,802,8,886,1045,0,100,62,25,6057,18,2446,2,0,2441,2216,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,1140,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:  1140<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU","count:     0<br />local_res_atom_non_h_electron_sum: 29<br />res_name: CU"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(133,173,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"CU","legendgroup":"CU","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9291,31249,55127,5825,3234,1244,12746,1620,1,8,3,870,802,8,886,1045,0,100,62,25,6057,18,2446,2,0,2441,2216,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,2,3362,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:  3362<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: DMS"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(114,176,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"DMS","legendgroup":"DMS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9285,16725,55127,5825,3234,1244,12746,1620,1,8,3,870,802,8,886,1045,0,100,62,25,6057,18,2446,2,0,2441,2216,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,6,14524,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     6<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count: 14524<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EDO"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(91,179,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"EDO","legendgroup":"EDO","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9285,16725,55106,5815,3217,1221,12732,740,1,8,3,870,802,8,886,1045,0,100,62,25,6057,18,2446,2,0,2441,2216,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,21,10,17,23,14,880,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:    21<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:    10<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:    17<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:    23<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:    14<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:   880<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: EPE"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(57,182,0,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"EPE","legendgroup":"EPE","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9285,16725,55106,5815,3217,1221,12732,740,1,8,3,870,802,8,886,1045,0,98,41,25,6057,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,21,0,0,2,0,0,0,7,2216,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:    21<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     7<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:  2216<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FAD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,184,31,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FAD","legendgroup":"FAD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9285,15991,55106,5815,3217,1221,12732,740,1,8,3,870,802,8,886,1045,0,98,41,25,6057,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,734,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:   734<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,186,66,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FE","legendgroup":"FE","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9285,15297,55106,5815,3217,1221,12732,740,1,8,3,870,802,8,886,1045,0,98,41,25,6057,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,694,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:   694<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2","count:     0<br />local_res_atom_non_h_electron_sum: 26<br />res_name: FE2"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,188,89,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FE2","legendgroup":"FE2","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9285,15297,55106,5815,3217,1221,12732,740,0,8,3,868,802,8,883,1,0,98,41,25,6057,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,0,0,1,0,0,2,0,0,3,1044,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:  1044<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMN"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,190,108,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FMN","legendgroup":"FMN","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9284,13970,55106,5815,3217,1221,12732,740,0,8,3,868,802,8,883,1,0,98,41,25,6057,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,1,1327,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:  1327<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: FMT"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,191,125,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"FMT","legendgroup":"FMT","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9284,13970,55106,5815,3217,1221,12731,740,0,6,3,867,802,8,89,1,0,98,41,25,6057,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,1,0,0,2,0,1,0,0,794,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:   794<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GDP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,192,141,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"GDP","legendgroup":"GDP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9283,13944,35970,5815,3217,1221,12731,740,0,6,3,867,802,8,89,1,0,98,41,25,6057,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,1,26,19136,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:    26<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count: 19136<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: GOL"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,193,156,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"GOL","legendgroup":"GOL","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9283,13944,35970,5815,3217,1221,12731,740,0,6,3,867,802,8,89,1,0,90,36,25,5323,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,8,5,0,734,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     8<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     5<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:   734<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEC"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,193,170,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"HEC","legendgroup":"HEC","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9283,13944,35970,5815,3217,1221,12731,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,2,25,5277,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:    25<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:  5277<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: HEM"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,192,184,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"HEM","legendgroup":"HEM","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9283,13944,35970,2441,3217,1221,12731,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,3374,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:  3374<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD","count:     0<br />local_res_atom_non_h_electron_sum: 53<br />res_name: IOD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,191,196,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"IOD","legendgroup":"IOD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,6966,13944,35970,2441,3217,1221,12731,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,2317,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:  2317<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K","count:     0<br />local_res_atom_non_h_electron_sum: 19<br />res_name: K"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,189,208,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"K","legendgroup":"K","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,6966,13944,35970,2439,2240,1094,12731,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,2,977,127,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:   977<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:   127<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MAN"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,187,219,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MAN","legendgroup":"MAN","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,6966,13944,35964,2433,2240,1091,11420,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,6,6,0,3,1311,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     6<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     6<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:  1311<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MES"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,184,229,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MES","legendgroup":"MES","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9,13944,35964,2433,2240,1091,11420,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,6957,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:  6957<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG","count:     0<br />local_res_atom_non_h_electron_sum: 12<br />res_name: MG"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,180,239,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MG","legendgroup":"MG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9,13923,35943,2210,837,1083,11420,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,21,21,223,1403,8,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:    21<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:    21<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:   223<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:  1403<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     8<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MLY"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,176,246,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MLY","legendgroup":"MLY","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9,11848,35943,2210,837,1083,11420,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,2075,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:  2075<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN","count:     0<br />local_res_atom_non_h_electron_sum: 25<br />res_name: MN"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,171,253,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MN","legendgroup":"MN","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9,11846,35939,836,837,1083,11420,740,0,6,3,867,802,8,89,0,0,90,34,0,46,16,2446,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,2,4,1374,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     4<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:  1374<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: MPD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(0,165,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"MPD","legendgroup":"MPD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1480,9,11846,35939,836,830,1083,11420,740,0,6,3,862,802,6,21,0,0,90,0,0,46,13,9,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,7,0,0,0,0,0,0,5,0,2,68,0,0,0,34,0,0,3,2437,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     7<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     5<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:    68<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:    34<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:  2437<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAD"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(82,158,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NAD","legendgroup":"NAD","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,9,11846,35939,835,829,1062,12,740,0,6,3,862,802,6,21,0,0,90,0,0,46,13,9,2,0,2434,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[1,0,0,0,1,1,21,11408,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:    21<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count: 11408<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAG"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(121,151,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NAG","legendgroup":"NAG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,9,11846,35939,835,829,1062,12,740,0,4,0,859,802,3,8,0,0,21,0,0,9,1,0,1,0,866,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,0,0,0,2,3,3,0,3,13,0,0,69,0,0,37,12,9,1,0,1568,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:    13<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:    69<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:    37<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:    12<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     9<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:  1568<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NAP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(149,144,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NAP","legendgroup":"NAP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,9,11846,35939,835,829,1062,12,740,0,4,0,859,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,0,0,0,0,0,0,3,3,8,0,0,21,0,0,9,1,0,1,0,866,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     8<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:    21<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     9<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:   866<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: NDP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(172,136,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NDP","legendgroup":"NDP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,9,10948,35939,835,829,1062,12,740,0,4,0,859,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,898,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:   898<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI","count:     0<br />local_res_atom_non_h_electron_sum: 28<br />res_name: NI"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(191,128,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NI","legendgroup":"NI","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,9,10068,35939,835,829,1062,12,740,0,4,0,859,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,880,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:   880<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3","count:     0<br />local_res_atom_non_h_electron_sum: 31<br />res_name: NO3"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(207,120,255,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"NO3","legendgroup":"NO3","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,7,10045,33743,835,829,1062,12,740,0,4,0,859,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,2,23,2196,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:    23<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:  2196<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PEG"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(220,113,250,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PEG","legendgroup":"PEG","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,7,10031,33641,801,715,1,12,740,0,4,0,859,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,14,102,34,114,1061,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:    14<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:   102<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:    34<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:   114<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:  1061<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PG4"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(231,107,243,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PG4","legendgroup":"PG4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,7,10014,33619,795,0,1,12,740,0,4,0,859,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,17,22,6,715,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:    17<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:    22<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     6<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:   715<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PGE"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(240,102,234,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PGE","legendgroup":"PGE","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,7,10014,33619,795,0,1,0,1,0,4,0,859,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,12,739,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:    12<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:   739<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PLP"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(247,99,224,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PLP","legendgroup":"PLP","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,5,10005,28343,795,0,1,0,1,0,4,0,859,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,2,9,5276,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     2<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     9<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:  5276<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: PO4"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(252,97,213,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"PO4","legendgroup":"PO4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,5,10005,28343,795,0,1,0,0,0,3,0,859,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,0,0,1,0,1,0,0,799,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:   799<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SAH"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,97,201,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"SAH","legendgroup":"SAH","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,5,10005,28343,795,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,0,0,0,0,1,0,0,0,3,0,859,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     1<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     3<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:   859<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SF4"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,98,188,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"SF4","legendgroup":"SF4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,0,9993,6,795,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,5,12,28337,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     5<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:    12<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count: 28337<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: SO4"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,101,174,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"SO4","legendgroup":"SO4","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[1479,0,9989,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,4,6,795,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     4<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     6<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:   795<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS","count:     0<br />local_res_atom_non_h_electron_sum: NA<br />res_name: TRS"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,104,159,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"TRS","legendgroup":"TRS","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[0,0,9989,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[1479,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:  1479<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX","count:     0<br />local_res_atom_non_h_electron_sum:  6<br />res_name: UNX"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(255,108,144,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"UNX","legendgroup":"UNX","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null},{"orientation":"v","width":[13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827587,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586,13.9310344827586],"base":[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"x":[0,13.9310344827586,27.8620689655172,41.7931034482759,55.7241379310345,69.6551724137931,83.5862068965517,97.5172413793103,111.448275862069,125.379310344828,139.310344827586,153.241379310345,167.172413793103,181.103448275862,195.034482758621,208.965517241379,222.896551724138,236.827586206897,250.758620689655,264.689655172414,278.620689655172,292.551724137931,306.48275862069,320.413793103448,334.344827586207,348.275862068966,362.206896551724,376.137931034483,390.068965517241,404],"y":[0,0,9989,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0],"text":["count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:  9989<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN","count:     0<br />local_res_atom_non_h_electron_sum: 30<br />res_name: ZN"],"type":"bar","marker":{"autocolorscale":false,"color":"rgba(252,113,127,1)","line":{"width":1.88976377952756,"color":"transparent"}},"name":"ZN","legendgroup":"ZN","showlegend":true,"xaxis":"x","yaxis":"y","hoverinfo":"text","frame":null}],"layout":{"margin":{"t":26.2283105022831,"r":7.30593607305936,"b":40.1826484018265,"l":54.7945205479452},"plot_bgcolor":"rgba(235,235,235,1)","paper_bgcolor":"rgba(255,255,255,1)","font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-27.8620689655172,431.862068965517],"tickmode":"array","ticktext":["0","100","200","300","400"],"tickvals":[0,100,200,300,400],"categoryorder":"array","categoryarray":["0","100","200","300","400"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"y","title":"local_res_atom_non_h_electron_sum","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"hoverformat":".2f"},"yaxis":{"domain":[0,1],"automargin":true,"type":"linear","autorange":false,"range":[-3067.55,64418.55],"tickmode":"array","ticktext":["0","20000","40000","60000"],"tickvals":[0,20000,40000,60000],"categoryorder":"array","categoryarray":["0","20000","40000","60000"],"nticks":null,"ticks":"outside","tickcolor":"rgba(51,51,51,1)","ticklen":3.65296803652968,"tickwidth":0.66417600664176,"showticklabels":true,"tickfont":{"color":"rgba(77,77,77,1)","family":"","size":11.689497716895},"tickangle":-0,"showline":false,"linecolor":null,"linewidth":0,"showgrid":true,"gridcolor":"rgba(255,255,255,1)","gridwidth":0.66417600664176,"zeroline":false,"anchor":"x","title":"count","titlefont":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"hoverformat":".2f"},"shapes":[{"type":"rect","fillcolor":null,"line":{"color":null,"width":0,"linetype":[]},"yref":"paper","xref":"paper","x0":0,"x1":1,"y0":0,"y1":1}],"showlegend":true,"legend":{"bgcolor":"rgba(255,255,255,1)","bordercolor":"transparent","borderwidth":1.88976377952756,"font":{"color":"rgba(0,0,0,1)","family":"","size":11.689497716895},"y":0.93503937007874},"annotations":[{"text":"res_name","x":1.02,"y":1,"showarrow":false,"ax":0,"ay":0,"font":{"color":"rgba(0,0,0,1)","family":"","size":14.6118721461187},"xref":"paper","yref":"paper","textangle":-0,"xanchor":"left","yanchor":"bottom","legendTitle":true}],"hovermode":"closest","barmode":"relative"},"config":{"doubleClick":"reset","modeBarButtonsToAdd":[{"name":"Collaborate","icon":{"width":1000,"ascent":500,"descent":-50,"path":"M487 375c7-10 9-23 5-36l-79-259c-3-12-11-23-22-31-11-8-22-12-35-12l-263 0c-15 0-29 5-43 15-13 10-23 23-28 37-5 13-5 25-1 37 0 0 0 3 1 7 1 5 1 8 1 11 0 2 0 4-1 6 0 3-1 5-1 6 1 2 2 4 3 6 1 2 2 4 4 6 2 3 4 5 5 7 5 7 9 16 13 26 4 10 7 19 9 26 0 2 0 5 0 9-1 4-1 6 0 8 0 2 2 5 4 8 3 3 5 5 5 7 4 6 8 15 12 26 4 11 7 19 7 26 1 1 0 4 0 9-1 4-1 7 0 8 1 2 3 5 6 8 4 4 6 6 6 7 4 5 8 13 13 24 4 11 7 20 7 28 1 1 0 4 0 7-1 3-1 6-1 7 0 2 1 4 3 6 1 1 3 4 5 6 2 3 3 5 5 6 1 2 3 5 4 9 2 3 3 7 5 10 1 3 2 6 4 10 2 4 4 7 6 9 2 3 4 5 7 7 3 2 7 3 11 3 3 0 8 0 13-1l0-1c7 2 12 2 14 2l218 0c14 0 25-5 32-16 8-10 10-23 6-37l-79-259c-7-22-13-37-20-43-7-7-19-10-37-10l-248 0c-5 0-9-2-11-5-2-3-2-7 0-12 4-13 18-20 41-20l264 0c5 0 10 2 16 5 5 3 8 6 10 11l85 282c2 5 2 10 2 17 7-3 13-7 17-13z m-304 0c-1-3-1-5 0-7 1-1 3-2 6-2l174 0c2 0 4 1 7 2 2 2 4 4 5 7l6 18c0 3 0 5-1 7-1 1-3 2-6 2l-173 0c-3 0-5-1-8-2-2-2-4-4-4-7z m-24-73c-1-3-1-5 0-7 2-2 3-2 6-2l174 0c2 0 5 0 7 2 3 2 4 4 5 7l6 18c1 2 0 5-1 6-1 2-3 3-5 3l-174 0c-3 0-5-1-7-3-3-1-4-4-5-6z"},"click":"function(gd) { \n        // is this being viewed in RStudio?\n        if (location.search == '?viewer_pane=1') {\n          alert('To learn about plotly for collaboration, visit:\\n https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html');\n        } else {\n          window.open('https://cpsievert.github.io/plotly_book/plot-ly-for-collaboration.html', '_blank');\n        }\n      }"}],"cloud":false},"source":"A","attrs":{"2348540d62ad":{"x":{},"fill":{},"type":"bar"}},"cur_data":"2348540d62ad","visdat":{"2348540d62ad":["function (y) ","x"]},"highlight":{"on":"plotly_click","persistent":false,"dynamic":false,"selectize":false,"opacityDim":0.2,"selected":{"opacity":1},"debounce":0},"base_url":"https://plot.ly"},"evals":["config.modeBarButtonsToAdd.0.click"],"jsHooks":[]}</script><!--/html_preserve-->

#Niezgodnoœæ liczby atomów

Tabela w tej sekcji przedstawia 10 klas z najwiêksz¹ niezgodnoœci¹ liczby atomów.


```r
niezg_atom <- data_final_50 %>%
  select(res_name, local_res_atom_non_h_count, dict_atom_non_h_count) %>%
  mutate(Roznica_atomow = dict_atom_non_h_count - local_res_atom_non_h_count) %>%
  group_by(res_name) %>%
  summarize(Niezgodnosc_atomow = sum(Roznica_atomow)) %>%
  arrange(desc(Niezgodnosc_atomow)) %>% head(10)

knitr::kable(niezg_atom) %>%
  kable_styling("striped", full_width = FALSE) %>% 
  column_spec(1, width = "2.5cm")
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> res_name </th>
   <th style="text-align:right;"> Niezgodnosc_atomow </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAG </td>
   <td style="text-align:right;"> 11195 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> CLA </td>
   <td style="text-align:right;"> 3052 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> 1PE </td>
   <td style="text-align:right;"> 2749 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MLY </td>
   <td style="text-align:right;"> 2343 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAP </td>
   <td style="text-align:right;"> 2182 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAD </td>
   <td style="text-align:right;"> 1833 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> COA </td>
   <td style="text-align:right;"> 1593 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> PG4 </td>
   <td style="text-align:right;"> 1239 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MAN </td>
   <td style="text-align:right;"> 986 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NDP </td>
   <td style="text-align:right;"> 763 </td>
  </tr>
</tbody>
</table>

#Niezgodnoœæ liczby elektronów

Tabela w tej sekcji przedstawia 10 klas z najwiêksz¹ niezgodnoœci¹ liczby elektronów.


```r
niezg_elect <- data_final_50 %>%
  select(res_name, local_res_atom_non_h_electron_sum, dict_atom_non_h_electron_sum) %>%
  mutate(Roznica_elektronow = dict_atom_non_h_electron_sum - local_res_atom_non_h_electron_sum) %>%
  group_by(res_name) %>%
  summarize(Niezgodnosc_elektronow = sum(Roznica_elektronow)) %>%
  arrange(desc(Niezgodnosc_elektronow)) %>% head(10)

knitr::kable(niezg_elect) %>%
  kable_styling("striped", full_width = FALSE) %>% 
  column_spec(1, width = "2.5cm")
```

<table class="table table-striped" style="width: auto !important; margin-left: auto; margin-right: auto;">
 <thead>
  <tr>
   <th style="text-align:left;"> res_name </th>
   <th style="text-align:right;"> Niezgodnosc_elektronow </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAG </td>
   <td style="text-align:right;"> 89528 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> 1PE </td>
   <td style="text-align:right;"> 18734 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> CLA </td>
   <td style="text-align:right;"> 18586 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MLY </td>
   <td style="text-align:right;"> 17441 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAP </td>
   <td style="text-align:right;"> 14817 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> NAD </td>
   <td style="text-align:right;"> 12302 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> COA </td>
   <td style="text-align:right;"> 11817 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> PG4 </td>
   <td style="text-align:right;"> 8372 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> MAN </td>
   <td style="text-align:right;"> 7888 </td>
  </tr>
  <tr>
   <td style="text-align:left;width: 2.5cm; "> PLP </td>
   <td style="text-align:right;"> 5232 </td>
  </tr>
</tbody>
</table>

#Wykresy part_01

Poni¿sze wykresy prezentuj¹ rozk³ad wartoœci kolumn 'part_01...' z podzia³em na klasy res_name.


```r
data_part <- data_final_50 %>% select(res_name,part_01_shape_segments_count:part_01_density_Z_4_0)
data_part_01 <- data_part %>% select(-res_name)
names_list <- as.list(colnames(data_part_01))

plot.fn.with.aes <- function(y.axis, x.axis) {
  
      ggplot(data_part, 
             aes_string(x = x.axis, 
                        y = y.axis,
                        color = x.axis)) + 
            geom_point() + theme(axis.text.x = element_blank()) +
            stat_summary(fun.y=mean, colour="black", geom="point", 
               shape=4, size=2, show.legend = FALSE) +
            stat_summary(fun.y=mean, colour="black", geom="text", show.legend = FALSE,
               hjust=-0.2, size=2.6, angle = 90, aes( label=round(..y.., digits=1)))
}

lapply(names_list, plot.fn.with.aes, x.axis = "res_name")
```

![](ZED_files/figure-html/part_01-1.png)<!-- -->![](ZED_files/figure-html/part_01-2.png)<!-- -->![](ZED_files/figure-html/part_01-3.png)<!-- -->![](ZED_files/figure-html/part_01-4.png)<!-- -->![](ZED_files/figure-html/part_01-5.png)<!-- -->![](ZED_files/figure-html/part_01-6.png)<!-- -->![](ZED_files/figure-html/part_01-7.png)<!-- -->![](ZED_files/figure-html/part_01-8.png)<!-- -->![](ZED_files/figure-html/part_01-9.png)<!-- -->![](ZED_files/figure-html/part_01-10.png)<!-- -->![](ZED_files/figure-html/part_01-11.png)<!-- -->![](ZED_files/figure-html/part_01-12.png)<!-- -->![](ZED_files/figure-html/part_01-13.png)<!-- -->![](ZED_files/figure-html/part_01-14.png)<!-- -->![](ZED_files/figure-html/part_01-15.png)<!-- -->![](ZED_files/figure-html/part_01-16.png)<!-- -->![](ZED_files/figure-html/part_01-17.png)<!-- -->![](ZED_files/figure-html/part_01-18.png)<!-- -->![](ZED_files/figure-html/part_01-19.png)<!-- -->![](ZED_files/figure-html/part_01-20.png)<!-- -->![](ZED_files/figure-html/part_01-21.png)<!-- -->![](ZED_files/figure-html/part_01-22.png)<!-- -->![](ZED_files/figure-html/part_01-23.png)<!-- -->![](ZED_files/figure-html/part_01-24.png)<!-- -->![](ZED_files/figure-html/part_01-25.png)<!-- -->![](ZED_files/figure-html/part_01-26.png)<!-- -->![](ZED_files/figure-html/part_01-27.png)<!-- -->![](ZED_files/figure-html/part_01-28.png)<!-- -->![](ZED_files/figure-html/part_01-29.png)<!-- -->![](ZED_files/figure-html/part_01-30.png)<!-- -->![](ZED_files/figure-html/part_01-31.png)<!-- -->![](ZED_files/figure-html/part_01-32.png)<!-- -->![](ZED_files/figure-html/part_01-33.png)<!-- -->![](ZED_files/figure-html/part_01-34.png)<!-- -->![](ZED_files/figure-html/part_01-35.png)<!-- -->![](ZED_files/figure-html/part_01-36.png)<!-- -->![](ZED_files/figure-html/part_01-37.png)<!-- -->![](ZED_files/figure-html/part_01-38.png)<!-- -->![](ZED_files/figure-html/part_01-39.png)<!-- -->![](ZED_files/figure-html/part_01-40.png)<!-- -->![](ZED_files/figure-html/part_01-41.png)<!-- -->![](ZED_files/figure-html/part_01-42.png)<!-- -->![](ZED_files/figure-html/part_01-43.png)<!-- -->![](ZED_files/figure-html/part_01-44.png)<!-- -->![](ZED_files/figure-html/part_01-45.png)<!-- -->![](ZED_files/figure-html/part_01-46.png)<!-- -->![](ZED_files/figure-html/part_01-47.png)<!-- -->![](ZED_files/figure-html/part_01-48.png)<!-- -->![](ZED_files/figure-html/part_01-49.png)<!-- -->![](ZED_files/figure-html/part_01-50.png)<!-- -->![](ZED_files/figure-html/part_01-51.png)<!-- -->![](ZED_files/figure-html/part_01-52.png)<!-- -->![](ZED_files/figure-html/part_01-53.png)<!-- -->![](ZED_files/figure-html/part_01-54.png)<!-- -->![](ZED_files/figure-html/part_01-55.png)<!-- -->![](ZED_files/figure-html/part_01-56.png)<!-- -->![](ZED_files/figure-html/part_01-57.png)<!-- -->![](ZED_files/figure-html/part_01-58.png)<!-- -->![](ZED_files/figure-html/part_01-59.png)<!-- -->![](ZED_files/figure-html/part_01-60.png)<!-- -->![](ZED_files/figure-html/part_01-61.png)<!-- -->![](ZED_files/figure-html/part_01-62.png)<!-- -->![](ZED_files/figure-html/part_01-63.png)<!-- -->![](ZED_files/figure-html/part_01-64.png)<!-- -->![](ZED_files/figure-html/part_01-65.png)<!-- -->![](ZED_files/figure-html/part_01-66.png)<!-- -->![](ZED_files/figure-html/part_01-67.png)<!-- -->![](ZED_files/figure-html/part_01-68.png)<!-- -->![](ZED_files/figure-html/part_01-69.png)<!-- -->![](ZED_files/figure-html/part_01-70.png)<!-- -->![](ZED_files/figure-html/part_01-71.png)<!-- -->![](ZED_files/figure-html/part_01-72.png)<!-- -->![](ZED_files/figure-html/part_01-73.png)<!-- -->![](ZED_files/figure-html/part_01-74.png)<!-- -->![](ZED_files/figure-html/part_01-75.png)<!-- -->![](ZED_files/figure-html/part_01-76.png)<!-- -->![](ZED_files/figure-html/part_01-77.png)<!-- -->![](ZED_files/figure-html/part_01-78.png)<!-- -->![](ZED_files/figure-html/part_01-79.png)<!-- -->![](ZED_files/figure-html/part_01-80.png)<!-- -->![](ZED_files/figure-html/part_01-81.png)<!-- -->![](ZED_files/figure-html/part_01-82.png)<!-- -->![](ZED_files/figure-html/part_01-83.png)<!-- -->![](ZED_files/figure-html/part_01-84.png)<!-- -->![](ZED_files/figure-html/part_01-85.png)<!-- -->![](ZED_files/figure-html/part_01-86.png)<!-- -->![](ZED_files/figure-html/part_01-87.png)<!-- -->![](ZED_files/figure-html/part_01-88.png)<!-- -->![](ZED_files/figure-html/part_01-89.png)<!-- -->![](ZED_files/figure-html/part_01-90.png)<!-- -->![](ZED_files/figure-html/part_01-91.png)<!-- -->![](ZED_files/figure-html/part_01-92.png)<!-- -->![](ZED_files/figure-html/part_01-93.png)<!-- -->![](ZED_files/figure-html/part_01-94.png)<!-- -->![](ZED_files/figure-html/part_01-95.png)<!-- -->![](ZED_files/figure-html/part_01-96.png)<!-- -->![](ZED_files/figure-html/part_01-97.png)<!-- -->![](ZED_files/figure-html/part_01-98.png)<!-- -->![](ZED_files/figure-html/part_01-99.png)<!-- -->![](ZED_files/figure-html/part_01-100.png)<!-- -->![](ZED_files/figure-html/part_01-101.png)<!-- -->![](ZED_files/figure-html/part_01-102.png)<!-- -->![](ZED_files/figure-html/part_01-103.png)<!-- -->![](ZED_files/figure-html/part_01-104.png)<!-- -->![](ZED_files/figure-html/part_01-105.png)<!-- -->![](ZED_files/figure-html/part_01-106.png)<!-- -->

#Regresja

Przeprowadzenie regresji oraz wyliczenie miar R^2 i RMSE dla przewidywania liczby atomów:


```r
idx <- createDataPartition(y = data_final_50$local_res_atom_non_h_count, p = .6, list = FALSE)
training <- data_final_50[idx, ]
testing <- data_final_50[-idx, ]
ctrl <- trainControl(method = "repeatedcv", number = 2,repeats = 5)

lm <- train(local_res_atom_non_h_count ~ .,
                data = training,
                method = "lm",
                metric = "RMSE",
                trControl = ctrl)
lmp<-predict(lm, newdata=testing)
postResample(lmp,testing$local_res_atom_non_h_count)
```

```
##       RMSE   Rsquared        MAE 
## 0.09388563 0.99994728 0.01515412
```

Przeprowadzenie regresji oraz wyliczenie miar R^2 i RMSE dla przewidywania liczby elektronów:


```r
idx <- createDataPartition(y = data_final_50$local_res_atom_non_h_electron_sum, p = .6, list = FALSE)
training <- data_final_50[idx, ]
testing <- data_final_50[-idx, ]
ctrl <- trainControl(method = "repeatedcv", number = 2,repeats = 5)

lm <- train(local_res_atom_non_h_electron_sum ~ .,
                data = training,
                method = "lm",
                metric = "RMSE",
                trControl = ctrl)
lmp<-predict(lm, newdata=testing)
postResample(lmp,testing$local_res_atom_non_h_electron_sum)
```

```
##      RMSE  Rsquared       MAE 
## 1.6640206 0.9996398 0.1031926
```

#Klasyfikator

Tworzenie klasyfikatora przewiduj¹cego wartoœæ atrybutu res_name:


```r
idx <- createDataPartition(data_final_50$res_name, p = .6, list = FALSE)

training <- data_final_50[idx,]

testing <- data_final_50[-idx,]

ctrl <- trainControl(method="repeatedcv", number=2, repeats = 5)

classifier <- train(as.factor(res_name) ~ .,
             data = training,
             method = "rf",
             trControl = ctrl,
             ntree = 4)

predictions <- predict(classifier, newdata = testing)

results <- confusionMatrix(data = predictions, 
                factor(testing[,1]))

knitr::kable(results$overall)
```

<table>
 <thead>
  <tr>
   <th style="text-align:left;">   </th>
   <th style="text-align:right;"> x </th>
  </tr>
 </thead>
<tbody>
  <tr>
   <td style="text-align:left;"> Accuracy </td>
   <td style="text-align:right;"> 0.9641502 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> Kappa </td>
   <td style="text-align:right;"> 0.9616800 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AccuracyLower </td>
   <td style="text-align:right;"> 0.9627634 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AccuracyUpper </td>
   <td style="text-align:right;"> 0.9654993 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AccuracyNull </td>
   <td style="text-align:right;"> 0.1581994 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> AccuracyPValue </td>
   <td style="text-align:right;"> 0.0000000 </td>
  </tr>
  <tr>
   <td style="text-align:left;"> McnemarPValue </td>
   <td style="text-align:right;"> NaN </td>
  </tr>
</tbody>
</table>
