#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 03:34:29 2017

@author: keidaiiiyama
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#乱数で生徒の志望学科リストを作成する
class Student:
    import numpy as np
    import random
    import matplotlib.pyplot as plt

    faculty_num = 7 #学科の数
    f_id = [0,1,2,3,4,5,6] #学科名
    def __init__(self,sn,seed):
        self.seed_num = seed
        self.np.random.seed(seed) #生徒の得点分布
        self.random.seed(seed) #生徒の学科希望
        self.student_num = sn * 100 #学生の数
        self.faculty_cap = [5*sn,8*sn,10*sn,12*sn,17*sn,22*sn,26*sn]#学科の定員
        self.min_score = [0,0,0,0,0,0,0] #第一段階最低点
        self.student_score_raw = self.np.random.normal(50,10,self.student_num)
        self.student_score = self.np.sort(self.student_score_raw) #生徒を得点順に並べたもの
        self.student_apply = [0,0,0,0,0,0,0] * self.student_num #生徒の学科希望順
        self.apply_position = [0] * self.student_num #その生徒が今何番目の志望先に申請しているか

    #乱数で生徒の志望学科リストを作成する
    def apply_list(self):
        for i in range(self.student_num):
            self.random.shuffle(self.f_id) #志望学科をランダムに設定
            self.student_apply[i] = [self.f_id[0],self.f_id[1],self.f_id[2],\
                               self.f_id[3],self.f_id[4],self.f_id[5],self.f_id[6]]

    def show_hist(self):
        fig1 = self.plt.figure()
        ax1 = fig1.add_subplot(1,1,1)
        ax1.hist(self.student_score, bins=50)
        ax1.set_title('student score histogram $\mu=50,\ \sigma=10$')
        ax1.set_xlabel('score')
        ax1.set_ylabel('number of students')
        fig1.show()

#第一段階まで
class FirstStage(Student):
    def __init__(self,sn,seed,ratio):
        Student.__init__(self,sn,seed)
        self.assigned = [[],[],[],[],[],[],[]] #各学科ごとの内定した学生(第一段階)
        self.not_assigned = [] #第一段階未内定者
        self.first_stage_ratio = ratio

    #定員の7割までを第一希望者内の成績で獲得
    def first_stage(self):
        for i in range(self.student_num):
            self.insert_student(i)

    def insert_student(self,student_id):
        first_apply = self.student_apply[student_id][self.apply_position[student_id]] #その生徒のapply_positon志望先
        insert_place = self.binary_search(self.assigned[first_apply],self.student_score[student_id])
        #点数足りず受け入れ拒否
        if(self.student_score[student_id] < self.min_score[first_apply]):
            self.not_assigned.append(student_id) #未内定者リストに追加
            return 0
        #点数足りる
        else:
            #学科リストに人を追加
            self.assigned[first_apply].insert(insert_place,student_id)
            if(student_id in self.not_assigned):
                self.not_assigned.remove(student_id) #未内定者リストから削除
            #追加により定員に達した場合、これ以降進学最低点が定義される
        if(len(self.assigned[first_apply]) == int(self.faculty_cap[first_apply] * self.first_stage_ratio)):
            if(len(self.assigned[first_apply]) == 0):
                return 0
            else:
                min_score_student = self.assigned[first_apply][0]
                self.min_score[first_apply] = self.student_score[min_score_student]
        #追加により定員を超えた場合、最低得点者は弾き出される
        elif(len(self.assigned[first_apply]) > int(self.faculty_cap[first_apply] * self.first_stage_ratio)):
            min_score_student = self.assigned[first_apply][0]
            self.not_assigned.append(min_score_student) #得点最下位者を未内定リストに追加
            self.assigned[first_apply].pop(0) #学科内定者名簿から最低得点者を削除
            min_score_student = self.assigned[first_apply][0] #進学最低点を更新
            self.min_score[first_apply] = self.student_score[min_score_student]

    #2分探索により、進学者名簿（得点低い順)に生徒を挿入する場所を計算する
    def binary_search(self,list,target):
        low = 0
        high = len(list)
        if high==0:
            return 0;
        if target < self.student_score[list[0]]:
            return 0
        else:
            while (high-low != 1):
                center = int((low+high)/2)
                if(self.student_score[list[center]] == target):
                    high = center + 1
                    break
                elif(self.student_score[list[center]] < target):
                    low = center
                elif(target < self.student_score[list[center]]):
                    high = center
        return high #targetはhigh番目（list[high-1])とhigh+1番目(list[high])の間にある

#第二段階まで
class SecondStage(FirstStage):
    def __init__(self,sn,seed,ratio,is_new):
        FirstStage.__init__(self,sn,seed,ratio)
        self.assigned2 = [[],[],[],[],[],[],[]] #第二段階内定者（全体)
        self.not_assigned2 = [] #第二段階未内定者(新規)
        self.not_assigned2_past = [] #前の判定での未内定者
        self.assigned2_stage = [[],[],[],[],[],[],[]] #第n志望順目における内定者
        self.first_assigned_rate = [0,0,0,0,0,0,0]
        self.is_new = is_new

    #第一段階の結果を受けた初期化処理
    def init_cap_nassigned(self):
                #faculty_capの更新
        for i in range(self.faculty_num):
            self.min_score[i] = 0
            self.first_assigned_rate[i] = len(self.assigned[i])/self.faculty_cap[i] * 100
            self.faculty_cap[i] = self.faculty_cap[i] - len(self.assigned[i])

        #not_assigned2_pastを生成
        for student in self.not_assigned:
            self.not_assigned2_past.append(student)

    #ソート用比較関数
    def key_func(self,n):
        return self.student_score[n]
    #各段階ごとの学科別進学最低点、最高点を表示
    def show_result(self):
        print("\n")
        print("結果")
        print("注:進学最低、最高点の後の数字は、その点数をとった人の番号を表す")
        for j in range(self.faculty_num):
            print("\n")
            print("学科番号:",j)
            #print("第1段階学科内定者：",assigned[j])
            print("第1段階内定者数",len(self.assigned[j]),"(",'%02.0f' % self.first_assigned_rate[j],"% )")
            if(self.first_assigned_rate[j] < self.first_stage_ratio *100):
                print("底割れ")
            print("第1段階進学最低点",self.student_score[self.assigned[j][0]],"(",self.assigned[j][0],")")
            print("第1段階進学最高点",self.student_score[self.assigned[j][len(self.assigned[j])-1]],
                "(",self.assigned[j][len(self.assigned[j])-1],")")
            #print("第2段階学科内定者：",assigned2[j])
            print("第2段階内定者数",len(self.assigned2[j]))
            min_score_student = min(self.assigned2[j],key = self.key_func)
            min_score = self.student_score[min_score_student]
            max_score_student = max(self.assigned2[j],key = self.key_func)
            max_score = self.student_score[max_score_student]
            print("第2段階進学最低点",min_score,"(",min_score_student,")")
            print("第2段階進学最高点",max_score,"(",max_score_student,")")
    #分布図を表示
    def show_figure(self):
        for i in range(len(self.apply_position)):
            self.apply_position[i] = self.apply_position[i] + 1

        #得点の順位で3つの集団に当分割
        x1 = self.student_score[0:int(self.student_num/5)]
        y1 = self.apply_position[0:int(self.student_num/5)]
        x2 = self.student_score[int(self.student_num/5):int(2*self.student_num/5)]
        y2 = self.apply_position[int(self.student_num/5):int(2*self.student_num/5)]
        x3 = self.student_score[int(2*self.student_num/5):int(3*self.student_num/5)]
        y3 = self.apply_position[int(2*self.student_num/5):int(3*self.student_num/5)]
        x4 = self.student_score[int(3*self.student_num/5):int(4*self.student_num/5)]
        y4 = self.apply_position[int(3*self.student_num/5):int(4*self.student_num/5)]
        x5 = self.student_score[int(4*self.student_num/5):int(self.student_num)]
        y5 = self.apply_position[int(4*self.student_num/5):int(self.student_num)]

        fig2 = self.plt.figure()
        ax = fig2.add_subplot(1,1,1)

        ax.scatter(x1,y1,s=2,c='red',label='bottom1/5')
        ax.scatter(x2,y2,s=2,c='blue',label='bottom2/5')
        ax.scatter(x3,y3,s=2,c='green',label='top3/5')
        ax.scatter(x4,y4,s=2,c='yellow',label='top2/5')
        ax.scatter(x5,y5,s=2,c='purple',label='top1/5')

        if(self.is_new==0):
            ax.set_title('score-order_of_choice scatter plot(old_system)')
            ax.set_xlabel('score')
            ax.set_ylabel('the choice order of the assigned faculty')
            fig2.show()
            filename = "shinfuri_old.png"
            self.plt.savefig(filename)
        if(self.is_new==1):
            ax.set_title('score-order_of_choice scatter plot(new system)')
            ax.set_xlabel('score')
            ax.set_ylabel('the choice order of the assigned faculty')
            fig2.show()
            filename = "shinfuri_new.png"
            self.plt.savefig(filename)

#旧進振りシステム
class OldShinfuri(SecondStage):
    def __init__(self,sn,seed,ratio):
        SecondStage.__init__(self,sn,seed,ratio,is_new=0)

    def simulate(self):
        self.apply_list()
        self.first_stage() #第一段階(共通)
        self.second_stage_old()
        self.show_figure()

    def second_stage_old(self):
        self.init_cap_nassigned()
        print("\n\n\n")
        print("乱数番号:",self.seed_num)
        print("第二段階(旧)アルゴリズム")
        order = 0
        #全員が学科に内定するまでループを回す
        while len(self.not_assigned2_past) > 0:
            print("第",order,"希望")
            for j in range(self.faculty_num):
                print("faculty",j,"min_score",self.min_score[j],"assigned",len(self.assigned2_stage[j]))
            print("残り定員",self.faculty_cap)
            print("not_assigned",len(self.not_assigned2_past))
            #格納用リストを空ける
            del self.not_assigned2[:]
            for k in range(self.faculty_num):
                del self.assigned2_stage[k][:]
            #内定名簿を更新
            for student in self.not_assigned2_past:
                self.insert_student2_old(student)

            del self.not_assigned2_past[:]    #更新のため、全削除
            for student in self.not_assigned2:
                self.apply_position[student] = self.apply_position[student] + 1 #内定しなかった人は希望を一つ下げる
                self.not_assigned2_past.append(student) #not_assigned_pastを更新
            for l in range(self.faculty_num):
                for student in self.assigned2_stage[l]:
                    self.assigned2[l].append(student) #assigned2_stageの学生をassigned2に追加
            for m in range(self.faculty_num):
                self.faculty_cap[m] = self.faculty_cap[m] - len(self.assigned2_stage[m]) #学科の残り定員を更新
            order = order + 1
        self.show_result()  #結果表示

    def insert_student2_old(self,student_id):
        #print("apply_position",apply_position[student_id])
        apply = self.student_apply[student_id][self.apply_position[student_id]] #その生徒のapply_positon志望先
        insert_place = self.binary_search(self.assigned2_stage[apply],self.student_score[student_id])
        #点数足りないかすでに前stageまでに定員超過で受け入れ拒否
        if(self.student_score[student_id] < self.min_score[apply] or self.faculty_cap[apply] == 0):
            self.not_assigned2.append(student_id) #未内定者リストに追加
        #点数足りる
        else:
            #学科リストに人を追加
            self.assigned2_stage[apply].insert(insert_place,student_id)
            #追加により定員に達した場合、これ以降進学最低点が定義される
            if(len(self.assigned2_stage[apply]) == self.faculty_cap[apply]):
                min_score_student = self.assigned2_stage[apply][0]
                self.min_score[apply] = self.student_score[min_score_student]
            #追加により定員を超えた場合、最低得点者は弾き出される
            elif(len(self.assigned2_stage[apply]) > self.faculty_cap[apply]):
                min_score_student = self.assigned2_stage[apply][0]
                self.not_assigned2.append(min_score_student) #得点最下位者を未内定リストに追加
                self.assigned2_stage[apply].pop(0) #stageでの学科内定者名簿から最低得点者を削除
                min_score_student = self.assigned2_stage[apply][0] #進学最低点を更新
                self.min_score[apply] = self.student_score[min_score_student]

#新進振りシステム
class NewShinfuri(SecondStage):
    def __init__(self,sn,seed,ratio):
        SecondStage.__init__(self,sn,seed,ratio,is_new=1)

    def simulate(self):
        self.apply_list()
        self.first_stage() #第一段階(共通)
        self.second_stage_new()
        self.show_figure()

    def second_stage_new(self):
        self.init_cap_nassigned()
        print("\n")
        print("第二段階(新)アルゴリズム")
        loop = 0
        #全員が学科に内定するまでループを回す
        while len(self.not_assigned2_past) > 0:
            print("第",loop,"ループ")
            for j in range(self.faculty_num):
                print("faculty",j,"min_score",self.min_score[j],"仮内定人数",len(self.assigned2[j]))
            print("not_assigned",len(self.not_assigned2_past))

            del self.not_assigned2[:]
            for student in self.not_assigned2_past:
                self.insert_student2_new(student) #内定名簿を更新
            del self.not_assigned2_past[:] #更新のため、全削除
            for student in self.not_assigned2:
                self.not_assigned2_past.append(student) #not_assigned_pastを更新
            loop = loop + 1
        self.show_result()  #結果表示

    def insert_student2_new(self,student_id):
        apply = self.student_apply[student_id][self.apply_position[student_id]] #その生徒のapply_position番目志望先
        insert_place = self.binary_search(self.assigned2[apply],self.student_score[student_id])
        #点数足りず受け入れ拒否
        if(self.student_score[student_id] < self.min_score[apply]):
            self.not_assigned2.append(student_id) #未内定者リストに追加
            self.apply_position[student_id] = self.apply_position[student_id] + 1 #内定しなかったので次は志望順位の1つ低いところへ
        #点数足りる
        else:
            #学科リストに人を追加
            self.assigned2[apply].insert(insert_place,student_id)
            #追加により定員に達した場合、これ以降進学最低点が定義される
            if(len(self.assigned2[apply]) == self.faculty_cap[apply]):
                min_score_student = self.assigned2[apply][0] #最低得点者
                self.min_score[apply] = self.student_score[min_score_student] #最低得点
            #追加により定員を超えた場合、最低得点者は弾き出される
            elif(len(self.assigned2[apply]) > self.faculty_cap[apply]):
                min_score_student = self.assigned2[apply][0] #得点最下位者(リストの先頭)
                self.not_assigned2.append(min_score_student) #得点最下位者を未内定リストに追加
                self.assigned2[apply].pop(0) #学科内定者名簿から最低得点者を削除
                min_score_student_new = self.assigned2[apply][0] #進学最低点を更新
                self.min_score[apply] = self.student_score[min_score_student_new]
                #apply_positionを更新
                self.apply_position[min_score_student] = self.apply_position[min_score_student] + 1

#引数：(生徒数、乱数の番号、第一段階で獲得する割合)
#比較するために新旧の進振りのシミュレーションの生徒作成には同じ番号の乱数を使う
class ShinfuriSimulator():
    def __init__(self,sn,sd,rt):
        self.number_of_students = sn
        self.seed = sd
        self.first_stage_ratio = rt

    def shinfuri_simulate(self):
        self.number_of_students = int(self.number_of_students/100)
        s_old = OldShinfuri(self.number_of_students,self.seed,self.first_stage_ratio)
        s_old.show_hist()
        s_old.simulate() #旧システム
        s_new = NewShinfuri(self.number_of_students,self.seed,self.first_stage_ratio)
        s_new.simulate() #新システム

if __name__ == "__main__":
    #乱数を変えて実験してみる
    ss = ShinfuriSimulator(3000,57,0.7)
    #ss2 = ShinfuriSimulator(3000,31,0.7)
    #ss3 = ShinfuriSimulator(3000,57,0.7)

    ss.shinfuri_simulate()
    #ss2.shinfuri_simulate()
    #ss3.shinfuri_simulate()
