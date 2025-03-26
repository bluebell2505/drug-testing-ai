
import React from 'react';
import Navbar from '@/components/Navbar';
import LoginForm from '@/components/LoginForm';
import Footer from '@/components/Footer';

const Index: React.FC = () => {
  return (
    <div className="flex flex-col min-h-screen bg-gradient-to-b from-gray-800 to-black">
      <div className="flex flex-1">
        {/* Left section - Visible only on md screens and above */}
        <div className="hidden lg:flex flex-col justify-center p-12 flex-1 relative">
          {/* Background image with reduced opacity */}
          <div className="absolute inset-0 z-0">
            <img
              src="/drug-discovery-image.png"
              alt="Drug Discovery Illustration"
              className="w-full h-full object-cover opacity-90"
            />
            <div className="absolute inset-0 bg-gradient-to-r from-black/90 to-transparent"></div>
          </div>
          
          <div className="relative z-10 max-w-lg">
            <h1 className="mb-6">
              <span className="flex text-5xl font-bold mb-2">
                <span className="text-blue-400">Drug</span>
                <span className="text-white">Discov</span>
              </span>
              <span className="block text-4xl font-semibold text-blue-300">The AI powered Drug discovery</span>
            </h1>
            <p className="text-xl text-gray-300 mb-8">
              Accelerate drug discovery with AI-driven molecule generation
            </p>
          </div>
        </div>
        
        {/* Right section - Login form */}
        <div className="flex flex-col justify-center items-center p-6 md:p-0 flex-1 bg-gradient-to-b from-gray-800 to-black">
          <div className="mb-8 block lg:hidden">
            <span className="flex text-4xl font-bold">
              <span className="text-blue-400">Drug</span>
              <span className="text-white">Discov</span>
            </span>
          </div>
          <LoginForm />
        </div>
      </div>
      
      <Footer />
    </div>
  );
};

export default Index;
