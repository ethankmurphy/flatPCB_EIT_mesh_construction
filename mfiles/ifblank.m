function filename = ifblank(filename)

%-------------------------------------------------------------------------------
% Set the decimal point character
for i = 1:length(filename)
   if strcmp(' ',filename(i))==1
       filename(i) = '_';
   end
end