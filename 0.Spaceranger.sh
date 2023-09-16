spaceranger-1.3.0/spaceranger count --id=Control_1 \
                   --transcriptome=refdata-gex-mm10-2020-A  \
                   --fastqs=Rawdata/Con_1 \
                   --sample=Con_1-1,Con_1-2 \
                   --image=Con_1/Con_1.jpg \
                   --slide=V11A19-112 \
                   --area=A1 \
                   --loupe-alignment=0.json/V11A19-112-A1.json

spaceranger-1.3.0/spaceranger count --id=Control_2 \
                   --transcriptome=refdata-gex-mm10-2020-A  \
                   --fastqs=Con_2 \
                   --sample=Con_2-1,Con_2-2 \
                   --image=Con_2/Con_2.jpg \
                   --slide=V11A19-112 \
                   --area=B1 \
                   --loupe-alignment=0.json/V11A19-112-B1.json

spaceranger-1.3.0/spaceranger count --id=AAN_2W_1 \
                   --transcriptome=refdata-gex-mm10-2020-A  \
                   --fastqs=AA_2w_1 \
                   --sample=AA_2w_1-1,AA_2w_1-2 \
                   --image=AA_2w_1/AA_2w_1.jpg \
                   --slide=V11A19-112 \
                   --area=C1 \
                   --loupe-alignment=0.json/V11A19-112-C1.json

spaceranger-1.3.0/spaceranger count --id=AAN_2W_2 \
                   --transcriptome=refdata-gex-mm10-2020-A  \
                   --fastqs=AA_2w_2 \
                   --sample=AA_2w_2-1,AA_2w_2-2 \
                   --image=AA_2w_2/AA_2w_2.jpg \
                   --slide=V11A19-112 \
                   --area=D1 \
                   --loupe-alignment=0.json/V11A19-112-D1.json

spaceranger-1.3.0/spaceranger count --id=AAN_4W_1 \
                   --transcriptome=refdata-gex-mm10-2020-A  \
                   --fastqs=Rawdata/AA_4W_C59/ \
                   --sample=AA_4W_C59-1,AA_4W_C59-2 \
                   --image=Rawdata/AA_4W_C59/AA_4W_C59.jpg \
                   --slide=V11Y03-339 \
                   --area=D1 \
                   --loupe-alignment=0.json/V11Y03-339-D1.json

spaceranger-1.3.0/spaceranger count --id=AAN_4W_2 \
                   --transcriptome=refdata-gex-mm10-2020-A  \
                   --fastqs=Rawdata/AA_4W_C60/ \
                   --sample=AA_4W_C60-1 \
                   --image=Rawdata/AA_4W_C60/AA_4W_C60.jpg \
                   --slide=V11Y03-338 \
                   --area=D1 \
                   --loupe-alignment=0.json/V11Y03-338-D1.json
