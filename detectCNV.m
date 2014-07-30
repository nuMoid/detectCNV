clc; clear all;clear;

%Add functions path to matlab search path
functionname='detectCNV.m'; functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
addpath([functiondir '/data'])
%% Parameters
% No. indivisual
num_ind = 12877;
% No.chromosome
num_chr = 1;
% threshold of derivative for smoothed data 
thresh = 0.04;
% minimal length of a variation
thresh_len = 25;

%% Loading and normalizing data 
%load
ind = num_ind;
num_ind = num2str(num_ind);
filename_start = ['C:\Users\Droid\Desktop\test\' num_ind '\start.mat'];
start = load(filename_start);start = [start.start];%start(1) = [];
filename_count = ['C:\Users\Droid\Desktop\test\' num_ind '\count.mat'];
count = load(filename_count);count = [count.count];%count(1) = [];
filename_ch = ['C:\Users\Droid\Desktop\test\' num_ind '\ch.mat'];
ch = load(filename_ch);ch = [ch.ch];
%
index = find(ch == num_chr);
count = count(index);start = start(index);
count_mean = 5000;
count_multi = count/count_mean;
% conversing noise value to avarage
thresh_noise = 4;
idx_noise = find(count_multi>thresh_noise);
count_multi(idx_noise) = 1;

%% Filtering and smoothing data
%count_smooth = mslowess(start,count_multi,'SPAN',15);
count_filt = medfilt2(count_multi,[5,1]);
count_smooth = mslowess(start,count_filt,'SPAN',10);
%subplot(1,2,1);plot(start,count_multi,'.');subplot(1,2,2);plot(start,count_smooth,'.');
count_diff = diff([0;count_smooth]);

figure; hold on
plot(start,count_multi,'.');
line(start,count_smooth,'Color','r');
line(start,count_diff,'Color','g');
% plot chromosome ideogram 
hs_cytobands = cytobandread('hs_cytoBand.txt')
set(gcf,'color','w')
chromosomeplot(hs_cytobands, num_chr, 'AddToPlot', gca,'unit', 1)
hold off

%% Detecting anormaly
idx = find(abs(count_diff)>thresh);

thresh_idx = 100;
j = 1;
for i = 1:length(idx)-1
    if idx(i+1)-idx(i)>thresh_idx
        idx_anor(j) = idx(i+1);
        j = j+1;
    end
end

%% Visualizing
figure; hold on
plot(start,count_multi,'.');

length_var = zeros(length(idx_anor),1);
mean = zeros(length(idx_anor),1);
for i = 1:length(idx_anor)-1
    head = idx_anor(i);
    tmp = find(idx == (idx_anor(i+1)));
    tail = idx(tmp-1);
    length_var(i) = tail-head;
    mean(i) = sum(count_multi(head:tail))/(tail-head+1);
    line(start(head:tail),mean(i), 'color','r');
    rectangle('Position',[start(head),-1,start(tail)-start(head),4]);
end

head = idx_anor(end);
tail = idx(end);
length_var(end) = tail-head;
mean(end) = sum(count_multi(head:tail))/(tail-head+1);
line(start(head:tail),mean(end), 'color','r');
rectangle('Position',[start(head),-1,start(tail)-start(head),4]);
% xlim([10^7,2*10^7]);
ylim([0,3]);
xlabel('Location');
ylabel('Read count');
hold off
%% Modifying
% remove false positive
j = 1;
for i = 1:length(idx_anor)
    if length_var(i)<12
        idx_s(j) = idx_anor(i);
        j = j+1;
    end
end

j = 1;
for i = 1:length(idx_anor)
    if length_var(i)<thresh_len
        idx_ts(j) = idx_anor(i);
        j = j+1;
    end
end

j = 1;
for i = 1:length(idx_anor)
    if mean(i)>0.8 && mean(i)<1.2
        idx_tf(j) = idx_anor(i);
        j = j+1;
    end
end

idx_fp = union(intersect(idx_tf,idx_ts),idx_s);

% remove region which is too close to centromere
% find centromere
idx_acen = strcmpi(hs_cytobands.GieStains, 'acen');
ends_acen = hs_cytobands.BandEndBPs(idx_acen);
position_acen  = ends_acen(1:2:end);
position = position_acen(num_chr)/10000;
length_rej = 1000;

j = 1;
for i = 1:length(idx_anor)
    if start(idx_anor(i))/10000>(position-length_rej) &&...
            start(idx_anor(i))/10000<(position+length_rej)
        idx_cen(j) = idx_anor(i);
        j = j+1;
    end
end

% remove region on both end
j = 1;
for i = 1:length(idx_anor)
    if start(idx_anor(i))/10000<1000 || start(idx_anor(i))/10000>length(index)-1000
        idx_end(j) = idx_anor(i);
        j = j+1;
    end
end

idx_tmp = union(idx_cen,idx_end);
idx_del = union(idx_tmp,idx_fp);
idx_anor2 = setdiff(idx_anor,idx_del);

%% Visualizing results after modifying
figure; hold on
plot(start,count_multi,'.');

length_var2 = zeros(length(idx_anor2),1);
for i = 1:length(idx_anor2)
    head = idx_anor2(i);
    tmp = find(idx_anor == head);
    tmp2 = find(idx == (idx_anor(tmp+1)));
    tail = idx(tmp2-1);
    length_var2(i) = tail-head;
    rectangle('Position',[start(head),-1,start(tail)-start(head),4]);
end

% head = idx_anor2(end);
% tmp = find(idx_anor == head);
% tmp2 = find(idx == (idx_anor(tmp+1)));
% tail = idx(tmp2-1);
% length_var(end) = tail-head;
% rectangle('Position',[start(head),-1,start(tail)-start(head),4]);
% xlim([10^7,2*10^7]);
ylim([0,3]);
xlabel('Location');
ylabel('Read count');
hold off
%% Outputting 
fprintf('%d copy number variations are found on Chr. %d of individual %d\n',...
    length(length_var2),num_chr,ind);