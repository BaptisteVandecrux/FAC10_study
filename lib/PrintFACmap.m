function [] = PrintFACmap(filename, map, mask, R)
    map(mask == 0) = NaN;
    geotiffwrite(filename, map,R,'CoordRefSysCode',3413)  
end