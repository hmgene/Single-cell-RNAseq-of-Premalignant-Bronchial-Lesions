for f in *.{r,sh}; do
    # 1. Edit the file and save the output to a temp file
    sed -e "s/good_repair/high_selfrenewal/g" \
        -e "s/bad_repair/low_selfrenewal/g" \
        -e "s/repairment/selfrenewal/g" "$f" > "$f.tmp"
        
    # 2. Overwrite the original file with the temporary one
    #mv "$f.tmp" "$f"
done
