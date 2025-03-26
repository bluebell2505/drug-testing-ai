
import React from 'react';
import Navbar from '@/components/Navbar';
import Footer from '@/components/Footer';

const Signup: React.FC = () => {
  return (
    <div className="min-h-screen bg-gray-900 flex flex-col">
      <Navbar />
      
      <div className="flex-1 flex items-center justify-center p-4">
        <h1 className="text-white text-2xl font-bold text-center">Coming Soon</h1>
      </div>
      
      <Footer />
    </div>
  );
};

export default Signup;
