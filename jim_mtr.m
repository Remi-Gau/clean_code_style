subjects = 2;

for ii =1:subjects
disp(['Processing subject' num2str(ii)]);
eval(['cd(''S' num2str(ii) ''');']);
jim_analyze;
FolderName = 'figures';
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  saveas(FigHandle, ['S' num2str(ii) '_' num2str(iFig) '.png']);
end
cd ..
end

disp('DONE');

close all; clear all;