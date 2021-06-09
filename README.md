# Matlab_Project_Jamming

  ## 🍏 monopulse radar 설계 및 jamming 모델링
  
  * monopulse radar
     * Sum-Difference Amplitude 방식으로 설계
     * beam을 4개의 part로 나누어서 각각의 beam에 들어오는 amplitude를 측정
     * 합차를 이용하여 목표물의 위치 계산


  * Jamming
     * Cross-Polarization 방식으로 각도기만 설계
     * 실제 목표물의 phase로부터 180[deg]만큼 차이가 나도록 함
     * phase를 azimuth angle과 elevation angle로 나누어 모델링  
    
---------------------------------------------


 ## 🍅 시뮬레이션 결과 분석
 
  1. monopulse radar modeling
    * 임의로 설정한 target의 azimuth angle과 elevation angle 
    * from 0 to -5[deg], -0.5[deg/s] 변화율
     
![4](https://user-images.githubusercontent.com/47622991/121314840-b9f81580-c942-11eb-89f3-8c45cb731174.PNG)
    
  2. 임의로 설정한 target angle 값과 추정한 target의 azi 및 elv의 angle 값
    * jamming을 걸지 않았을 경우
    * AWGN channel (mu = 0 sigma = 0.043) 환경에서 적용

![5](https://user-images.githubusercontent.com/47622991/121315992-d34d9180-c943-11eb-8efd-48d00eb7a76c.PNG)

  3. RMSE 측정 (no jamming)
    * no jamming 에서의 RMSE 측정 (실제 target의 위치와 얼마나 차이가 나는지)
    * azimuth angle 측정 (elevation angle 측정은 생략)

![6](https://user-images.githubusercontent.com/47622991/121316380-3c350980-c944-11eb-961a-c9db23f2b823.PNG)

  4. jamming vs no jamming
    * jamming을 걸어준 경우와 jamming이 없는 경우를 동시에 비교
    * cross-polarization 기법을 적용하여 실제 target의 angle 대비 180 [deg] 만큼 차이 발생

![7](https://user-images.githubusercontent.com/47622991/121316785-a3eb5480-c944-11eb-9522-51fe7528612b.PNG)

  5. RMSE 측정 (jamming)
    * jamming 상태에서의 RMSE 측정 (실제 target의 위치와 얼마나 차이가 나는지)
    * azimuth angle 측정 (elevation angle 측정은 생략)

![8](https://user-images.githubusercontent.com/47622991/121317179-02183780-c945-11eb-9e7e-53ced0675ee5.PNG)

  6. RMSE 측정 (jamming의 세기)
    * no jamming vs jamming (A=1) vs jamming (A=3)
    * jamming의 세기에 따라 RMSE 측정이 어떻게 변화 하는지 비교

![9](https://user-images.githubusercontent.com/47622991/121317409-3a1f7a80-c945-11eb-81a1-92200c8507ac.PNG)


