//Code specific for Arduino Mega 2560
// PORTE 4 is line 2
// PORTH 0 is line 17

float duty_cycle_off = 0.5;
int line_period_us = 60;
int off_time = line_period_us * duty_cycle_off;
int a = off_time/3;
int c = a*2;

void setup() {
  // put your setup code here, to run once:
DDRE = B11111100;        //Pin 1 of PORTD is an input, all others are outputs
DDRH = B00000000;        //Pin 1 of PORTD is an input, all others are outputs


DDRH = B00000000; 
PORTH = B11111111; 

noInterrupts();

}

void loop() {
  
if ((PINH & (B00000001))==1){
  PORTE = (1<<PD4);    //Pin 2 of portd as now the logic value 1
  delayMicroseconds(a); //16
  PORTE = (0<<PD4);    //Pin 2 of portd as now the logic value 0
  delayMicroseconds(off_time); //34
  PORTE = (1<<PD4);    //Pin 2 of portd as now the logic value 1
  delayMicroseconds(c); //32
  PORTE = (0<<PD4);    //Pin 2 of portd as now the logic value 0
  delayMicroseconds(off_time);  //34
  PORTE = (1<<PD4);    //Pin 2 of portd as now the logic value 1
}

}
