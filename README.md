# Matlab_Project_Jamming

  ## ๐ monopulse radar ์ค๊ณ ๋ฐ jamming ๋ชจ๋ธ๋ง
  
  * monopulse radar
     * Sum-Difference Amplitude ๋ฐฉ์์ผ๋ก ์ค๊ณ
     * beam์ 4๊ฐ์ part๋ก ๋๋์ด์ ๊ฐ๊ฐ์ beam์ ๋ค์ด์ค๋ amplitude๋ฅผ ์ธก์ 
     * ํฉ์ฐจ๋ฅผ ์ด์ฉํ์ฌ ๋ชฉํ๋ฌผ์ ์์น ๊ณ์ฐ


  * Jamming
     * Cross-Polarization ๋ฐฉ์์ผ๋ก ๊ฐ๋๊ธฐ๋ง ์ค๊ณ
     * ์ค์  ๋ชฉํ๋ฌผ์ phase๋ก๋ถํฐ 180[deg]๋งํผ ์ฐจ์ด๊ฐ ๋๋๋ก ํจ
     * phase๋ฅผ azimuth angle๊ณผ elevation angle๋ก ๋๋์ด ๋ชจ๋ธ๋ง  
    
---------------------------------------------


 ## ๐ ์๋ฎฌ๋ ์ด์ ๊ฒฐ๊ณผ ๋ถ์
 
  **1. monopulse radar modeling**
 
  * ์์๋ก ์ค์ ํ target์ azimuth angle๊ณผ elevation angle 
  * from 0 to -5[deg], -0.5[deg/s] ๋ณํ์จ
     
![4](https://user-images.githubusercontent.com/47622991/121314840-b9f81580-c942-11eb-89f3-8c45cb731174.PNG)
    
  **2. ์์๋ก ์ค์ ํ target angle ๊ฐ๊ณผ ์ถ์ ํ target์ azi ๋ฐ elv์ angle ๊ฐ**

  * jamming์ ๊ฑธ์ง ์์์ ๊ฒฝ์ฐ
  * AWGN channel (mu = 0 sigma = 0.043) ํ๊ฒฝ์์ ์ ์ฉ

![5](https://user-images.githubusercontent.com/47622991/121315992-d34d9180-c943-11eb-8efd-48d00eb7a76c.PNG)

  **3. RMSE ์ธก์  (no jamming)**

  * no jamming ์์์ RMSE ์ธก์  (์ค์  target์ ์์น์ ์ผ๋ง๋ ์ฐจ์ด๊ฐ ๋๋์ง)
  * azimuth angle ์ธก์  (elevation angle ์ธก์ ์ ์๋ต)

![6](https://user-images.githubusercontent.com/47622991/121316380-3c350980-c944-11eb-961a-c9db23f2b823.PNG)

  **4. jamming vs no jamming**
  
  * jamming์ ๊ฑธ์ด์ค ๊ฒฝ์ฐ์ jamming์ด ์๋ ๊ฒฝ์ฐ๋ฅผ ๋์์ ๋น๊ต
  * cross-polarization ๊ธฐ๋ฒ์ ์ ์ฉํ์ฌ ์ค์  target์ angle ๋๋น 180 [deg] ๋งํผ ์ฐจ์ด ๋ฐ์

![7](https://user-images.githubusercontent.com/47622991/121316785-a3eb5480-c944-11eb-9522-51fe7528612b.PNG)

  **5. RMSE ์ธก์  (jamming)**
  
  * jamming ์ํ์์์ RMSE ์ธก์  (์ค์  target์ ์์น์ ์ผ๋ง๋ ์ฐจ์ด๊ฐ ๋๋์ง)
  * azimuth angle ์ธก์  (elevation angle ์ธก์ ์ ์๋ต)

![8](https://user-images.githubusercontent.com/47622991/121317179-02183780-c945-11eb-9e7e-53ced0675ee5.PNG)

  **6. RMSE ์ธก์  (jamming์ ์ธ๊ธฐ)**

  * no jamming vs jamming (A=1) vs jamming (A=3)
  * jamming์ ์ธ๊ธฐ์ ๋ฐ๋ผ RMSE ์ธก์ ์ด ์ด๋ป๊ฒ ๋ณํ ํ๋์ง ๋น๊ต

![9](https://user-images.githubusercontent.com/47622991/121317409-3a1f7a80-c945-11eb-81a1-92200c8507ac.PNG)


