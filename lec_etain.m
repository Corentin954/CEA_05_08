clear , close all;



resu=load("pkb.txt");
[~,n]=size(resu);

v=resu(:,1);
pkb=resu(:,2);


figure
hold on
plot(v,pkb)
xlabel('vol sp�cifique');
ylabel('Pkb');
grid();
hold off