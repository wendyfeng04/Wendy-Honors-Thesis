Hi! you can follow the instructions below to regenerate all my results:) 
# Step 1: Processing climate model inputs ($\lambda$, $\alpha$, and $E_0$)from FLUXNET data
(1) To process raw FLUXNET data, please refer to [Fluxdata_Pro.R](./Fluxdata_Pro.R), the input file should contain half-hourly FLUXNET data (FLUXNET_FLUXMET_HH).

(2) Once all raw FLUXNET data have been processed, run [Growing_Season_Sum.R](./Growing_Season_Sum.R) to generate a summary table of the growing-season parameters ($\lambda$, $\alpha$, and $E_0$), which will serve as model inputs. The script also calculates the average growing-season FLUXNET evapotranspiration ($ET$) for reference.

(3) A summary table of all FLUXNET-derived climatic inputs, along with other site-specific model parameters and information, including citations to their sources, can be found in [Sites_Data.csv](./Sites_Data.csv).

# Step 2: Plotting
To generate all the figures, please run the [Plot.ipynb](./Plot.ipynb) notebook in the same folder/environment with [sswm.py](./sswm.py), [growth.py](./growth.py), [param.py](./param.py),and [Sites_Data.csv](./Sites_Data.csv)

Please note that a substantial amount of the code was adapted from https://github.com/maoyab/OWUS 
