function S = txtToStruct(textfile)

  fid = fopen(textfile,'r');
  textline = fgets(fid);
  headers = strsplit(strtrim(textline),{' ','\t'});
  data = fscanf(fid, '%f',[length(headers) Inf]);
  fclose(fid);
  for hh = 1:length(headers)
    S.(headers{hh}) = data(hh,:);
  end

end