//Code specific for Arduino Mega 2560
// PORTE 4 is line 2
// PORTH 0 is line 17

int potPin = 2; // select the analog input pin for the pot
int val; // var to store incoming signal from sensor

int a = 10230; // 1023 * 10
int b = 38874; // 1023 * 38
int c = 30690; // 1023 * 30
int d = 38874; // 1023 * 38
int aa;
int bb;
int cc;
int dd;

void setup() {
  // put your setup code here, to run once:
DDRE = B11111100;        //Pin 1 of PORTD is an input, all others are outputs
DDRH = B00000000;        //Pin 1 of PORTD is an input, all others are outputs


DDRH = B00000000; 
PORTH = B11111111; 

val = analogRead(A0); // between 0 and 1023
aa = a / val;
bb = b / val;
cc = c / val;
dd = d / val;

noInterrupts();

}

void loop() {
  
if ((PINH & (B00000001))==1){
  
  
  PORTE = (1<<PD4);    //Pin 2 of portd as now the logic value 1
  delayMicroseconds(aa); //16
  PORTE = (0<<PD4);    //Pin 2 of portd as now the logic value 0
  delayMicroseconds(bb); //34
  PORTE = (1<<PD4);    //Pin 2 of portd as now the logic value 1
  delayMicroseconds(cc); //32
  PORTE = (0<<PD4);    //Pin 2 of portd as now the logic value 0
  delayMicroseconds(dd);  //34
  PORTE = (1<<PD4);    //Pin 2 of portd as now the logic value 1
}

}
